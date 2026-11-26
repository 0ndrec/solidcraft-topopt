//! A small, safe layer over the CUDA driver API, loaded at run time (`nvcuda.dll` /
//! `libcuda.so.1`): no CUDA toolkit is needed to build or run, only an NVIDIA driver to use a GPU.

use std::ffi::{CStr, CString, c_char, c_int, c_uint, c_void};
use std::sync::Arc;

type CuResult = c_int;
type CuDevice = c_int;
type CuContext = *mut c_void;
type CuModule = *mut c_void;
type CuFunction = *mut c_void;
type CuDevicePtr = u64;

const ATTR_MULTIPROCESSOR_COUNT: c_int = 16;
const ATTR_COMPUTE_CAPABILITY_MAJOR: c_int = 75;
const ATTR_COMPUTE_CAPABILITY_MINOR: c_int = 76;
const JIT_ERROR_LOG_BUFFER: c_int = 5;
const JIT_ERROR_LOG_BUFFER_SIZE_BYTES: c_int = 6;

/// The driver library, opened once and never unloaded.
#[cfg(windows)]
mod os {
    use std::ffi::{CStr, c_char, c_void};
    pub const LIBRARY: &CStr = c"nvcuda.dll";
    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn LoadLibraryA(name: *const c_char) -> *mut c_void;
        fn GetProcAddress(module: *mut c_void, name: *const c_char) -> *mut c_void;
    }
    pub fn open() -> *mut c_void {
        unsafe { LoadLibraryA(LIBRARY.as_ptr()) }
    }
    pub fn symbol(lib: *mut c_void, name: &CStr) -> *mut c_void {
        unsafe { GetProcAddress(lib, name.as_ptr()) }
    }
}

#[cfg(unix)]
mod os {
    use std::ffi::{CStr, c_char, c_int, c_void};
    pub const LIBRARY: &CStr = c"libcuda.so.1";
    const RTLD_NOW: c_int = 2;
    #[link(name = "dl")]
    unsafe extern "C" {
        fn dlopen(name: *const c_char, flags: c_int) -> *mut c_void;
        fn dlsym(lib: *mut c_void, name: *const c_char) -> *mut c_void;
    }
    pub fn open() -> *mut c_void {
        unsafe { dlopen(LIBRARY.as_ptr(), RTLD_NOW) }
    }
    pub fn symbol(lib: *mut c_void, name: &CStr) -> *mut c_void {
        unsafe { dlsym(lib, name.as_ptr()) }
    }
}

macro_rules! api {
    ($($name:ident: fn($($arg:ty),*);)*) => {
        #[allow(non_snake_case)]
        struct Api {
            $($name: unsafe extern "system" fn($($arg),*) -> CuResult,)*
        }

        impl Api {
            fn load() -> Result<Api, String> {
                let library = os::LIBRARY.to_string_lossy();
                let lib = os::open();
                if lib.is_null() {
                    return Err(format!("no NVIDIA driver ({library} not found)"));
                }
                Ok(Api {
                    $($name: {
                        let name = CString::new(stringify!($name)).unwrap();
                        let f = os::symbol(lib, &name);
                        if f.is_null() {
                            return Err(format!("{library} has no {} (driver too old?)", stringify!($name)));
                        }
                        // SAFETY: the driver exports this function with this signature.
                        unsafe { std::mem::transmute::<*mut c_void, unsafe extern "system" fn($($arg),*) -> CuResult>(f) }
                    },)*
                })
            }
        }
    };
}

api! {
    cuInit: fn(c_uint);
    cuDriverGetVersion: fn(*mut c_int);
    cuGetErrorString: fn(CuResult, *mut *const c_char);
    cuDeviceGetCount: fn(*mut c_int);
    cuDeviceGet: fn(*mut CuDevice, c_int);
    cuDeviceGetName: fn(*mut c_char, c_int, CuDevice);
    cuDeviceTotalMem_v2: fn(*mut usize, CuDevice);
    cuDeviceGetAttribute: fn(*mut c_int, c_int, CuDevice);
    cuDevicePrimaryCtxRetain: fn(*mut CuContext, CuDevice);
    cuDevicePrimaryCtxRelease_v2: fn(CuDevice);
    cuCtxSetCurrent: fn(CuContext);
    cuCtxSynchronize: fn();
    cuMemGetInfo_v2: fn(*mut usize, *mut usize);
    cuModuleLoadDataEx: fn(*mut CuModule, *const c_void, c_uint, *mut c_int, *mut *mut c_void);
    cuModuleUnload: fn(CuModule);
    cuModuleGetFunction: fn(*mut CuFunction, CuModule, *const c_char);
    cuMemAlloc_v2: fn(*mut CuDevicePtr, usize);
    cuMemFree_v2: fn(CuDevicePtr);
    cuMemcpyHtoD_v2: fn(CuDevicePtr, *const c_void, usize);
    cuMemcpyDtoH_v2: fn(*mut c_void, CuDevicePtr, usize);
    cuMemcpyDtoD_v2: fn(CuDevicePtr, CuDevicePtr, usize);
    cuMemsetD8_v2: fn(CuDevicePtr, u8, usize);
    cuFuncGetAttribute: fn(*mut c_int, c_int, CuFunction);
    cuLaunchKernel: fn(CuFunction, c_uint, c_uint, c_uint, c_uint, c_uint, c_uint, c_uint, *mut c_void, *mut *mut c_void, *mut *mut c_void);
}

impl Api {
    fn check(&self, code: CuResult, what: &str) -> Result<(), String> {
        if code == 0 {
            return Ok(());
        }
        let mut s: *const c_char = std::ptr::null();
        let text = if unsafe { (self.cuGetErrorString)(code, &mut s) } == 0 && !s.is_null() {
            unsafe { CStr::from_ptr(s) }.to_string_lossy().into_owned()
        } else {
            "unknown error".into()
        };
        Err(format!("CUDA {what}: {text} ({code})"))
    }
}

/// The driver, loaded and initialized.
pub struct Driver {
    api: Api,
    /// e.g. 13000 for CUDA 13.0.
    pub version: i32,
}

/// A GPU as the driver reports it.
#[derive(Clone, Debug, PartialEq)]
pub struct DeviceInfo {
    pub index: usize,
    pub name: String,
    pub memory: usize,
    pub compute_capability: (i32, i32),
    pub multiprocessors: i32,
}

impl Driver {
    pub fn load() -> Result<Arc<Driver>, String> {
        let api = Api::load()?;
        api.check(unsafe { (api.cuInit)(0) }, "init")?;
        let mut version = 0;
        api.check(unsafe { (api.cuDriverGetVersion)(&mut version) }, "driver version")?;
        Ok(Arc::new(Driver { api, version }))
    }

    pub fn devices(&self) -> Result<Vec<DeviceInfo>, String> {
        let a = &self.api;
        let mut count = 0;
        a.check(unsafe { (a.cuDeviceGetCount)(&mut count) }, "device count")?;
        (0..count.max(0) as usize).map(|i| self.device_info(i)).collect()
    }

    fn device_info(&self, index: usize) -> Result<DeviceInfo, String> {
        let a = &self.api;
        let mut dev = 0;
        a.check(unsafe { (a.cuDeviceGet)(&mut dev, index as c_int) }, "device")?;
        let mut name = [0 as c_char; 256];
        a.check(unsafe { (a.cuDeviceGetName)(name.as_mut_ptr(), 255, dev) }, "device name")?;
        let mut memory = 0usize;
        a.check(unsafe { (a.cuDeviceTotalMem_v2)(&mut memory, dev) }, "device memory")?;
        let attr = |which| -> Result<i32, String> {
            let mut v = 0;
            a.check(unsafe { (a.cuDeviceGetAttribute)(&mut v, which, dev) }, "device attribute")?;
            Ok(v)
        };
        Ok(DeviceInfo {
            index,
            name: unsafe { CStr::from_ptr(name.as_ptr()) }.to_string_lossy().trim().to_string(),
            memory,
            compute_capability: (attr(ATTR_COMPUTE_CAPABILITY_MAJOR)?, attr(ATTR_COMPUTE_CAPABILITY_MINOR)?),
            multiprocessors: attr(ATTR_MULTIPROCESSOR_COUNT)?,
        })
    }
}

/// A context on one GPU with the kernels loaded.
pub struct Gpu {
    driver: Arc<Driver>,
    device: CuDevice,
    ctx: CuContext,
    module: CuModule,
    functions: Vec<CuFunction>,
    pub info: DeviceInfo,
}

// SAFETY: the driver API is thread-safe; every entry point makes the context current first.
unsafe impl Send for Gpu {}
unsafe impl Sync for Gpu {}

/// One kernel argument.
#[derive(Clone, Copy, Debug)]
pub enum Arg {
    Ptr(u64),
    U32(u32),
    F32(f32),
    F64(f64),
    /// A parameter passed by value (a struct or an array): its bytes, which must stay valid
    /// until [`Gpu::launch`] returns (the driver copies them at launch).
    Bytes(*const c_void),
}

impl Gpu {
    /// Opens device `index` and loads `ptx`, resolving the kernels `names` (indexed as given).
    pub fn open(driver: Arc<Driver>, index: usize, ptx: &str, names: &[&str]) -> Result<Arc<Gpu>, String> {
        let info = driver.device_info(index)?;
        let a = &driver.api;
        let mut device = 0;
        a.check(unsafe { (a.cuDeviceGet)(&mut device, index as c_int) }, "device")?;
        let mut ctx = std::ptr::null_mut();
        a.check(unsafe { (a.cuDevicePrimaryCtxRetain)(&mut ctx, device) }, "context")?;
        a.check(unsafe { (a.cuCtxSetCurrent)(ctx) }, "context")?;
        let mut gpu = Gpu { driver: driver.clone(), device, ctx, module: std::ptr::null_mut(), functions: Vec::new(), info };
        let src = CString::new(ptx).map_err(|_| "the PTX contains a NUL byte".to_string())?;
        let mut log = vec![0u8; 8192];
        let mut opts = [JIT_ERROR_LOG_BUFFER, JIT_ERROR_LOG_BUFFER_SIZE_BYTES];
        let mut vals = [log.as_mut_ptr() as *mut c_void, log.len() as *mut c_void];
        let code = unsafe { (a.cuModuleLoadDataEx)(&mut gpu.module, src.as_ptr() as *const c_void, 2, opts.as_mut_ptr(), vals.as_mut_ptr()) };
        if code != 0 {
            let text = String::from_utf8_lossy(&log[..log.iter().position(|&b| b == 0).unwrap_or(0)]).into_owned();
            return Err(a.check(code, "loading the kernels").unwrap_err() + if text.is_empty() { "" } else { ": " } + &text);
        }
        for name in names {
            let mut f = std::ptr::null_mut();
            let c = CString::new(*name).unwrap();
            a.check(unsafe { (a.cuModuleGetFunction)(&mut f, gpu.module, c.as_ptr()) }, &format!("kernel {name}"))?;
            gpu.functions.push(f);
        }
        Ok(Arc::new(gpu))
    }

    fn api(&self) -> &Api {
        &self.driver.api
    }

    /// Makes this GPU's context current on the calling thread.
    pub fn bind(&self) -> Result<(), String> {
        self.api().check(unsafe { (self.api().cuCtxSetCurrent)(self.ctx) }, "context")
    }

    /// Free and total device memory, bytes.
    pub fn memory(&self) -> Result<(usize, usize), String> {
        self.bind()?;
        let (mut free, mut total) = (0, 0);
        self.api().check(unsafe { (self.api().cuMemGetInfo_v2)(&mut free, &mut total) }, "memory info")?;
        Ok((free, total))
    }

    pub fn synchronize(&self) -> Result<(), String> {
        self.api().check(unsafe { (self.api().cuCtxSynchronize)() }, "synchronize")
    }

    /// Registers per thread and local memory bytes per thread of kernel `k` (as the driver
    /// compiled it for this GPU).
    pub fn kernel_resources(&self, k: usize) -> Result<(i32, i32), String> {
        let a = self.api();
        let (mut regs, mut local) = (0, 0);
        a.check(unsafe { (a.cuFuncGetAttribute)(&mut regs, 4, self.functions[k]) }, "kernel attribute")?;
        a.check(unsafe { (a.cuFuncGetAttribute)(&mut local, 3, self.functions[k]) }, "kernel attribute")?;
        Ok((regs, local))
    }

    /// Launches kernel `k` on `blocks × threads` threads.
    pub fn launch(&self, k: usize, blocks: u32, threads: u32, args: &[Arg]) -> Result<(), String> {
        let mut vals: Vec<Arg> = args.to_vec();
        let mut ptrs: Vec<*mut c_void> = vals
            .iter_mut()
            .map(|a| match a {
                Arg::Ptr(v) => v as *mut u64 as *mut c_void,
                Arg::U32(v) => v as *mut u32 as *mut c_void,
                Arg::F32(v) => v as *mut f32 as *mut c_void,
                Arg::F64(v) => v as *mut f64 as *mut c_void,
                Arg::Bytes(p) => *p as *mut c_void,
            })
            .collect();
        let code = unsafe {
            (self.api().cuLaunchKernel)(self.functions[k], blocks.max(1), 1, 1, threads, 1, 1, 0, std::ptr::null_mut(), ptrs.as_mut_ptr(), std::ptr::null_mut())
        };
        self.api().check(code, "kernel launch")
    }
}

impl Drop for Gpu {
    fn drop(&mut self) {
        let a = &self.driver.api;
        unsafe {
            (a.cuCtxSetCurrent)(self.ctx);
            if !self.module.is_null() {
                (a.cuModuleUnload)(self.module);
            }
            (a.cuDevicePrimaryCtxRelease_v2)(self.device);
        }
    }
}

/// Device memory holding `len` values of `T`.
pub struct Buffer<T: Copy> {
    gpu: Arc<Gpu>,
    ptr: u64,
    len: usize,
    _t: std::marker::PhantomData<T>,
}

impl<T: Copy> Buffer<T> {
    pub fn new(gpu: &Arc<Gpu>, len: usize) -> Result<Self, String> {
        let a = gpu.api();
        let mut ptr = 0;
        let bytes = (len.max(1)) * std::mem::size_of::<T>();
        a.check(unsafe { (a.cuMemAlloc_v2)(&mut ptr, bytes) }, &format!("allocating {} MB", bytes.div_ceil(1 << 20)))?;
        let b = Buffer { gpu: gpu.clone(), ptr, len, _t: std::marker::PhantomData };
        b.zero()?;
        Ok(b)
    }

    pub fn from_slice(gpu: &Arc<Gpu>, data: &[T]) -> Result<Self, String> {
        let b = Self::new(gpu, data.len())?;
        b.upload(data)?;
        Ok(b)
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn arg(&self) -> Arg {
        Arg::Ptr(self.ptr)
    }

    /// The device address of element `i`.
    pub fn at(&self, i: usize) -> Arg {
        Arg::Ptr(self.ptr + (i * std::mem::size_of::<T>()) as u64)
    }

    fn bytes(&self) -> usize {
        self.len * std::mem::size_of::<T>()
    }

    pub fn upload(&self, data: &[T]) -> Result<(), String> {
        assert_eq!(data.len(), self.len, "upload size");
        let a = self.gpu.api();
        a.check(unsafe { (a.cuMemcpyHtoD_v2)(self.ptr, data.as_ptr() as *const c_void, self.bytes()) }, "upload")
    }

    pub fn download(&self, out: &mut [T]) -> Result<(), String> {
        assert_eq!(out.len(), self.len, "download size");
        let a = self.gpu.api();
        a.check(unsafe { (a.cuMemcpyDtoH_v2)(out.as_mut_ptr() as *mut c_void, self.ptr, self.bytes()) }, "download")
    }

    /// The first `out.len()` values.
    pub fn download_prefix(&self, out: &mut [T]) -> Result<(), String> {
        assert!(out.len() <= self.len, "download size");
        let a = self.gpu.api();
        let bytes = std::mem::size_of_val(out);
        a.check(unsafe { (a.cuMemcpyDtoH_v2)(out.as_mut_ptr() as *mut c_void, self.ptr, bytes) }, "download")
    }

    /// Value `i`.
    pub fn get(&self, i: usize) -> Result<T, String> {
        assert!(i < self.len, "index");
        let a = self.gpu.api();
        let mut v = std::mem::MaybeUninit::<T>::uninit();
        let at = self.ptr + (i * std::mem::size_of::<T>()) as u64;
        a.check(unsafe { (a.cuMemcpyDtoH_v2)(v.as_mut_ptr() as *mut c_void, at, std::mem::size_of::<T>()) }, "download")?;
        Ok(unsafe { v.assume_init() })
    }

    pub fn copy_from(&self, src: &Buffer<T>) -> Result<(), String> {
        assert_eq!(src.len, self.len, "copy size");
        let a = self.gpu.api();
        a.check(unsafe { (a.cuMemcpyDtoD_v2)(self.ptr, src.ptr, self.bytes()) }, "copy")
    }

    pub fn zero(&self) -> Result<(), String> {
        let a = self.gpu.api();
        a.check(unsafe { (a.cuMemsetD8_v2)(self.ptr, 0, self.bytes()) }, "memset")
    }
}

impl<T: Copy> Drop for Buffer<T> {
    fn drop(&mut self) {
        let a = self.gpu.api();
        unsafe {
            (a.cuCtxSetCurrent)(self.gpu.ctx);
            (a.cuMemFree_v2)(self.ptr);
        }
    }
}
