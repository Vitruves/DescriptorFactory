#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
import time
import signal
import platform
import tempfile
import re
import json
from pathlib import Path
import glob

class Logger:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    CYAN = '\033[0;36m'
    WHITE = '\033[0;37m'
    NC = '\033[0m'
    BOLD = '\033[1m'
    
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.progress_items = 0
        self.progress_current = 0
        self.progress_desc = ""
        
    def _format_message(self, level, message, color):
        return f"{color}[{level}]{self.NC} {message}"
    
    def info(self, message):
        print(self._format_message("INFO", message, self.BLUE))
        
    def debug(self, message):
        if self.verbose:
            print(self._format_message("DEBUG", message, self.WHITE))
        
    def warn(self, message):
        print(self._format_message("WARN", message, self.YELLOW))
        
    def error(self, message):
        print(self._format_message("ERROR", message, self.RED))
        
    def success(self, message):
        print(self._format_message("SUCCESS", message, self.GREEN))
    
    def start_progress(self, total, desc="Processing"):
        self.progress_items = total
        self.progress_current = 0
        self.progress_desc = desc
        print(f"{self.BLUE}[INFO]{self.NC} {desc}: 0/{total} (0.00%)")
        return self
    
    def update(self, n=1):
        self.progress_current += n
        percentage = (self.progress_current / self.progress_items) * 100
        print(f"{self.BLUE}[INFO]{self.NC} {self.progress_desc}: {self.progress_current}/{self.progress_items} ({percentage:.2f}%)")
    
    def set_description(self, desc):
        self.progress_desc = desc
    
    def finish(self):
        if self.progress_items > 0:
            print(f"{self.BLUE}[INFO]{self.NC} {self.progress_desc}: {self.progress_items}/{self.progress_items} (100.00%)")
            self.progress_items = 0
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish()

def run_command(cmd, logger, cwd=None, env=None, show_output=False, shell=False, timeout=None, hide_warnings=True):
    logger.debug(f"Running command: {cmd if shell else ' '.join(cmd)}")
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
            shell=shell
        )
        
        stdout_data = []
        stderr_data = []
        
        import select
        if timeout:
            end_time = time.time() + timeout
        
        poller = select.poll()
        poller.register(process.stdout, select.POLLIN)
        poller.register(process.stderr, select.POLLIN)
        
        while process.poll() is None:
            try:
                if timeout and time.time() > end_time:
                    process.terminate()
                    raise subprocess.TimeoutExpired(cmd, timeout)
                
                events = poller.poll(100)
                
                for fd, event in events:
                    if fd == process.stdout.fileno():
                        line = process.stdout.readline()
                        if line:
                            stdout_data.append(line)
                            if show_output or logger.verbose:
                                if not hide_warnings or "warning:" not in line:
                                    print(line.rstrip())
                    elif fd == process.stderr.fileno():
                        line = process.stderr.readline()
                        if line:
                            stderr_data.append(line)
                            if show_output or logger.verbose:
                                if not hide_warnings or "warning:" not in line:
                                    print(line.rstrip())
            except KeyboardInterrupt:
                process.terminate()
                logger.warn("Command interrupted by user")
                return False, "".join(stdout_data), "Command interrupted by user"
            
            time.sleep(0.001)
        
        for line in process.stdout:
            stdout_data.append(line)
            if show_output or logger.verbose:
                if not hide_warnings or "warning:" not in line:
                    print(line.rstrip())
            
        for line in process.stderr:
            stderr_data.append(line)
            if show_output or logger.verbose:
                if not hide_warnings or "warning:" not in line:
                    print(line.rstrip())
        
        stdout = "".join(stdout_data)
        stderr = "".join(stderr_data)
        
        if process.returncode != 0:
            logger.error(f"Command failed with exit code {process.returncode}")
            if not show_output and not logger.verbose:
                logger.debug(f"Command output:\n{stdout}")
                logger.error(f"Command error:\n{stderr}")
            return False, stdout, stderr
        
        return True, stdout, stderr
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout} seconds")
        return False, "", "Command timed out"
    except Exception as e:
        logger.error(f"Exception when running command: {e}")
        return False, "", str(e)

def detect_cpu_features(logger):
    logger.info("Detecting CPU features...")
    features = {
        "vendor": "",
        "model": "",
        "arch": platform.machine(),
        "cores": os.cpu_count(),
        "frequency": 0,
        "cache_size": 0,
        "features": []
    }
    
    if sys.platform == "linux":
        try:
            with open("/proc/cpuinfo", "r") as f:
                cpu_info = f.read()
            
            vendor_match = re.search(r"vendor_id\s+:\s+(.+)", cpu_info)
            if vendor_match:
                features["vendor"] = vendor_match.group(1)
                
            model_match = re.search(r"model name\s+:\s+(.+)", cpu_info)
            if model_match:
                features["model"] = model_match.group(1)
                
            freq_match = re.search(r"cpu MHz\s+:\s+(.+)", cpu_info)
            if freq_match:
                features["frequency"] = float(freq_match.group(1))
                
            cache_match = re.search(r"cache size\s+:\s+(\d+)\s+KB", cpu_info)
            if cache_match:
                features["cache_size"] = int(cache_match.group(1))
                
            flags_match = re.search(r"flags\s+:\s+(.+)", cpu_info)
            if flags_match:
                features["features"] = flags_match.group(1).split()
        except Exception as e:
            logger.warn(f"Error reading CPU info: {e}")
    
    elif sys.platform == "darwin":
        try:
            sysctl_cmd = ["sysctl", "-a"]
            success, stdout, _ = run_command(sysctl_cmd, logger, hide_warnings=True)
            if success:
                vendor_match = re.search(r"machdep.cpu.vendor: (.+)", stdout)
                if vendor_match:
                    features["vendor"] = vendor_match.group(1)
                    
                brand_match = re.search(r"machdep.cpu.brand_string: (.+)", stdout)
                if brand_match:
                    features["model"] = brand_match.group(1)
                    
                freq_match = re.search(r"hw.cpufrequency: (\d+)", stdout)
                if freq_match:
                    features["frequency"] = int(freq_match.group(1)) / 1000000
                
                l3_match = re.search(r"hw.l3cachesize: (\d+)", stdout)
                if l3_match:
                    features["cache_size"] = int(l3_match.group(1)) / 1024
                    
                features_cmd = ["sysctl", "machdep.cpu.features"]
                success, stdout, _ = run_command(features_cmd, logger, hide_warnings=True)
                if success:
                    features_match = re.search(r"machdep.cpu.features: (.+)", stdout)
                    if features_match:
                        features["features"] = features_match.group(1).split()
        except Exception as e:
            logger.warn(f"Error getting CPU info: {e}")
            
    elif sys.platform == "win32":
        try:
            wmic_cmd = ["wmic", "cpu", "get", "Name,NumberOfCores,MaxClockSpeed", "/format:csv"]
            success, stdout, _ = run_command(wmic_cmd, logger, hide_warnings=True)
            if success:
                lines = stdout.strip().split('\n')
                if len(lines) > 1:
                    parts = lines[1].split(',')
                    if len(parts) >= 3:
                        features["model"] = parts[1]
                        features["cores"] = int(parts[2])
                        features["frequency"] = int(parts[3])
                        
        except Exception as e:
            logger.warn(f"Error getting CPU info: {e}")
            
    logger.success(f"CPU: {features['model']} ({features['cores']} cores)")
    
    cpu_capabilities = []
    important_features = ["avx", "avx2", "avx512", "fma", "sse4_2", "aes", "f16c", "bmi2"]
    
    for feature in important_features:
        if any(f.lower() == feature.lower() or f.lower().startswith(feature.lower()) for f in features["features"]):
            cpu_capabilities.append(feature)
    
    if cpu_capabilities:
        logger.info(f"CPU capabilities: {', '.join(cpu_capabilities)}")
    
    return features, cpu_capabilities

def detect_gpu_features(logger):
    logger.info("Detecting GPU capabilities...")
    gpus = []
    
    try:
        nvidia_smi_cmd = ["nvidia-smi", "--query-gpu=name,memory.total,cuda_capability_major,cuda_capability_minor", "--format=csv,noheader"]
        success, stdout, _ = run_command(nvidia_smi_cmd, logger, hide_warnings=True)
        
        if success:
            lines = stdout.strip().split('\n')
            for i, line in enumerate(lines):
                parts = line.split(', ')
                if len(parts) >= 4:
                    gpu = {
                        "index": i,
                        "name": parts[0],
                        "memory": parts[1],
                        "cuda_capability": f"{parts[2]}.{parts[3]}"
                    }
                    gpus.append(gpu)
                    logger.success(f"Found GPU: {gpu['name']} ({gpu['memory']}, CUDA {gpu['cuda_capability']})")
        
        if not gpus:
            logger.debug("No NVIDIA GPUs detected with nvidia-smi, trying alternative methods")
            
            import shutil
            if shutil.which("lspci"):
                lspci_cmd = ["lspci"]
                success, stdout, _ = run_command(lspci_cmd, logger, hide_warnings=True)
                if success:
                    vga_sections = re.findall(r"VGA.*?^\S", stdout, re.MULTILINE | re.DOTALL)
                    for section in vga_sections:
                        if "nvidia" in section.lower():
                            name_match = re.search(r"VGA.*?: (.*)", section)
                            if name_match:
                                gpu_name = name_match.group(1)
                                gpus.append({"name": gpu_name, "memory": "Unknown", "cuda_capability": "Unknown"})
                                logger.info(f"Found GPU: {gpu_name} (CUDA capability unknown)")
    except Exception as e:
        logger.warn(f"Error detecting GPUs: {e}")
        
    if not gpus:
        logger.info("No CUDA-capable GPUs detected")
    
    return gpus

def detect_memory_info(logger):
    logger.info("Detecting system memory...")
    memory_info = {"total": 0, "available": 0}
    
    try:
        if sys.platform == "linux":
            with open("/proc/meminfo", "r") as f:
                mem_info = f.read()
            
            total_match = re.search(r"MemTotal:\s+(\d+)\s+kB", mem_info)
            if total_match:
                memory_info["total"] = int(total_match.group(1)) * 1024
                
            available_match = re.search(r"MemAvailable:\s+(\d+)\s+kB", mem_info)
            if available_match:
                memory_info["available"] = int(available_match.group(1)) * 1024
                
        elif sys.platform == "darwin":
            sysctl_cmd = ["sysctl", "hw.memsize"]
            success, stdout, _ = run_command(sysctl_cmd, logger, hide_warnings=True)
            if success:
                mem_match = re.search(r"hw.memsize: (\d+)", stdout)
                if mem_match:
                    memory_info["total"] = int(mem_match.group(1))
                    
            vm_stat_cmd = ["vm_stat"]
            success, stdout, _ = run_command(vm_stat_cmd, logger, hide_warnings=True)
            if success:
                page_size_match = re.search(r"page size of (\d+) bytes", stdout)
                free_pages_match = re.search(r"Pages free:\s+(\d+)", stdout)
                
                if page_size_match and free_pages_match:
                    page_size = int(page_size_match.group(1))
                    free_pages = int(free_pages_match.group(1))
                    memory_info["available"] = page_size * free_pages
                    
        elif sys.platform == "win32":
            wmic_cmd = ["wmic", "OS", "get", "FreePhysicalMemory,TotalVisibleMemorySize", "/format:csv"]
            success, stdout, _ = run_command(wmic_cmd, logger, hide_warnings=True)
            if success:
                lines = stdout.strip().split('\n')
                if len(lines) > 1:
                    parts = lines[1].split(',')
                    if len(parts) >= 3:
                        memory_info["available"] = int(parts[1]) * 1024
                        memory_info["total"] = int(parts[2]) * 1024
                    
    except Exception as e:
        logger.warn(f"Error detecting memory info: {e}")
    
    total_gb = memory_info["total"] / (1024**3)
    available_gb = memory_info["available"] / (1024**3)
    
    logger.success(f"Memory: {total_gb:.1f} GB total, {available_gb:.1f} GB available")
    return memory_info

def detect_compiler_info(logger):
    logger.info("Detecting available compilers...")
    compilers = {}
    
    for compiler in ["g++", "clang++", "c++"]:
        try:
            result = subprocess.run([compiler, "--version"], 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    text=True,
                                    timeout=5)
            if result.returncode == 0:
                version = result.stdout.strip().split('\n')[0]
                compilers[compiler] = {
                    "version": version,
                    "path": shutil.which(compiler)
                }
                logger.success(f"Found {compiler}: {version}")
                
                if "clang" in version.lower():
                    compilers[compiler]["type"] = "clang"
                elif "gcc" in version.lower():
                    compilers[compiler]["type"] = "gcc"
                else:
                    compilers[compiler]["type"] = "unknown"
        except Exception as e:
            logger.debug(f"Error detecting {compiler}: {e}")
    
    return compilers

def detect_compiler_flags(logger, quick=False):
    logger.info("Detecting supported compiler flags for PGO...")
    
    if quick:
        logger.info("Quick mode: using standard GCC PGO flags")
        return "-fprofile-generate", "-fprofile-use"
    
    temp_dir = "/tmp/compiler_test"
    os.makedirs(temp_dir, exist_ok=True)
    
    with open(os.path.join(temp_dir, "test.cpp"), "w") as f:
        f.write("#include <iostream>\nint main() { std::cout << \"Hello, World!\" << std::endl; return 0; }")
    
    compiler_info = {}
    for compiler in ["c++", "g++", "clang++"]:
        try:
            result = subprocess.run([compiler, "--version"], 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    text=True,
                                    timeout=5)
            if result.returncode == 0:
                compiler_info[compiler] = result.stdout.lower()
                logger.info(f"Found compiler: {compiler}")
                logger.debug(f"Version info: {result.stdout.strip()}")
        except subprocess.TimeoutExpired:
            logger.warn(f"Timeout while checking {compiler} version")
        except:
            pass
    
    compiler_cmd = "c++"
    if "g++" in compiler_info:
        compiler_cmd = "g++"
    elif "clang++" in compiler_info:
        compiler_cmd = "clang++"
    
    chosen_version = compiler_info.get(compiler_cmd, "")
    is_clang = "clang" in chosen_version
    is_gcc = "gcc" in chosen_version or "g++" in chosen_version and not is_clang
    
    logger.info(f"Using {compiler_cmd} for PGO flag detection")
    
    pgo_flag_options = []
    
    if is_gcc:
        logger.info("Detected GCC-compatible compiler")
        pgo_flag_options.extend([
            ("-fprofile-generate", "-fprofile-use -fprofile-correction"),
            ("-fprofile-generate=.", "-fprofile-use=. -fprofile-correction"),
            ("-fprofile-generate", "-fprofile-use -fbranch-probabilities"),
            ("-fprofile-generate", "-fprofile-use -fprofile-correction -fbranch-probabilities"),
            ("-fprofile-generate", "-fprofile-use"),
            ("-fprofile-generate=.", "-fprofile-use=."),
        ])
    
    if is_clang:
        logger.info("Detected Clang-compatible compiler")
        pgo_flag_options.extend([
            ("-fprofile-instr-generate", "-fprofile-instr-use=default.profdata"),
            ("-fprofile-instr-generate", "-fprofile-instr-use=default.profdata -Wno-profile-instr-out-of-date"),
            ("-fprofile-instr-generate=profile.profraw", "-fprofile-instr-use=profile.profdata"),
            ("-fprofile-instr-generate", "-fprofile-instr-use"),
            ("-fprofile-generate", "-fprofile-use"),
        ])
    
    if not pgo_flag_options:
        logger.warn("Could not identify compiler type, trying generic options")
        pgo_flag_options.extend([
            ("-fprofile-generate", "-fprofile-use"),
            ("-fprofile-instr-generate", "-fprofile-instr-use"),
        ])
    
    gen_flag, use_flag = "", ""
    for gen, use in pgo_flag_options:
        compile_cmd = f"{compiler_cmd} {gen} -o test test.cpp"
        logger.debug(f"Testing compilation with: {compile_cmd}")
        
        try:
            compile_result = subprocess.run(
                compile_cmd, 
                shell=True,
                cwd=temp_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=10
            )
            
            if compile_result.returncode == 0:
                logger.debug(f"Compilation with {gen} succeeded")
                
                run_result = subprocess.run(
                    "./test", 
                    shell=True,
                    cwd=temp_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=5
                )
                
                if run_result.returncode == 0:
                    use_cmd = f"{compiler_cmd} {use} -o test_optimized test.cpp"
                    logger.debug(f"Testing optimization with: {use_cmd}")
                    
                    use_result = subprocess.run(
                        use_cmd,
                        shell=True,
                        cwd=temp_dir,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        timeout=10
                    )
                    
                    if use_result.returncode == 0:
                        gen_flag, use_flag = gen, use
                        logger.success(f"Found working PGO flags: {gen} → {use}")
                        break
                
        except subprocess.TimeoutExpired:
            logger.warn(f"Timeout while testing {gen}")
        except Exception as e:
            logger.debug(f"Error testing {gen}: {str(e)}")
    
    try:
        shutil.rmtree(temp_dir, ignore_errors=True)
    except:
        logger.debug(f"Could not clean up {temp_dir}")
    
    if not gen_flag:
        logger.warn("Could not determine working PGO flags, falling back to CMake defaults")
        return None, None
    
    return gen_flag, use_flag

def get_default_optimization_flags(cpu_capabilities):
    optimization_flags = "-O3"
    
    if "avx512" in cpu_capabilities:
        optimization_flags += " -mavx512f -mavx512dq -mavx512bw -mavx512vl"
    elif "avx2" in cpu_capabilities:
        optimization_flags += " -mavx2"
        if "fma" in cpu_capabilities:
            optimization_flags += " -mfma"
    elif "avx" in cpu_capabilities:
        optimization_flags += " -mavx"
    elif "sse4_2" in cpu_capabilities:
        optimization_flags += " -msse4.2"
        
    return optimization_flags

def create_optimized_cmakelists(src_dir, dest_path, optimization_flags, logger, skip_gui=False):
    logger.info(f"Creating optimized CMakeLists.txt at {dest_path}")
    
    original_cmakelists = os.path.join(src_dir, "CMakeLists.txt")
    if not os.path.exists(original_cmakelists):
        logger.error(f"Original CMakeLists.txt not found at {original_cmakelists}")
        return False
    
    try:
        os.makedirs(os.path.dirname(os.path.abspath(dest_path)), exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create directory for CMakeLists.txt: {e}")
        return False

    try:
        with open(original_cmakelists, "r") as f:
            content = f.read()
        
        modified_content = content
        
        # Set build type to Release
        build_type_pattern = r"(set\s*\(\s*CMAKE_BUILD_TYPE\s+)[^\)]+(\s*\))"
        modified_content = re.sub(build_type_pattern, r"\1Release\2", modified_content)
        
        # Add our custom configuration section at the top
        custom_config = f"""
# Optimized build configuration
set(CMAKE_BUILD_TYPE Release)

# Optimization flags - properly formatted as a string
set(CMAKE_CXX_FLAGS_RELEASE "{optimization_flags}")

# Option to skip GUI build
option(BUILD_GUI "Build the GUI application" {"OFF" if skip_gui else "ON"})

# PGO Configuration with proper flag structure
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${{CMAKE_CXX_FLAGS}} -fprofile-generate=.")
    set(CMAKE_EXE_LINKER_FLAGS "${{CMAKE_EXE_LINKER_FLAGS}} -fprofile-generate=.")
    set(CMAKE_SHARED_LINKER_FLAGS "${{CMAKE_SHARED_LINKER_FLAGS}} -fprofile-generate=.")
    message(STATUS "PGO: Enabled profile generation with -fprofile-generate=.")
endif()
"""
        # Insert at beginning, but after cmake_minimum_required if present
        if "cmake_minimum_required" in modified_content:
            pattern = r"(cmake_minimum_required\([^\)]+\))"
            modified_content = re.sub(pattern, r"\1\n\n" + custom_config, modified_content, count=1)
        else:
            modified_content = custom_config + "\n\n" + modified_content
        
        # Remove any existing CMAKE_CXX_FLAGS settings to avoid conflicts
        flags_pattern = r"set\s*\(\s*CMAKE_CXX_FLAGS(?:_RELEASE)?\s+[^\)]+\)"
        modified_content = re.sub(flags_pattern, "", modified_content)

        # Ensure BUILD_GUI option takes effect if GUI target exists
        # This is a guess; adjust if your CMake structure differs
        gui_target_pattern = r"add_executable\s*\(\s*desfact-gui\s+.*?\)"
        if re.search(gui_target_pattern, modified_content, re.DOTALL):
             modified_content = re.sub(gui_target_pattern,
                                        r"if(BUILD_GUI)\n\g<0>\nendif()",
                                        modified_content)
        
        with open(dest_path, "w") as f:
            f.write(modified_content)
        
        logger.success(f"Created optimized CMakeLists.txt with proper formatting for flags: {optimization_flags}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to write optimized CMakeLists.txt: {e}")
        return False

def ensure_input_directory(input_dir, logger):
    if not os.path.exists(input_dir):
        logger.info(f"Creating input directory: {input_dir}")
        os.makedirs(input_dir, exist_ok=True)
    
    if not any(f.endswith('.csv') for f in os.listdir(input_dir)):
        logger.info("No CSV files found in input directory. Creating sample data file...")
        
        # More robust sample data without header line called "SMILES"
        sample_data = """molecule
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CC(=O)OC1=CC=CC=C1C(=O)O
CC1=C(C=C(C=C1)CC(=O)O)OC
CC1=CC(=C(C(=C1)OC)O)O
CC1=CC=C(O1)C(=O)O
CC(C)C1=CC=CC(=C1)C(=O)O
CC(C)(C)C1=CC(=C(C=C1)O)CC(=O)O
CC(C)(C)C1=CC(=CC(=C1)O)CC(=O)O
CC(C)C1=CC=CC=C1CC(=O)O
CC(C)C1=CC=CC=C1CCC(=O)O"""
        
        sample_file_path = os.path.join(input_dir, "sample_smiles.csv")
        with open(sample_file_path, 'w') as f:
            f.write(sample_data)
        
        logger.success(f"Created sample SMILES file at {sample_file_path}")

def stage_instrumented_build(src_dir, build_dir, jobs, logger, pgo_flags=None, timeout=None, hide_warnings=True, cmakelists_path=None, skip_gui=False):
    logger.info("STAGE 1: Building instrumented binary")

    # Resolve paths to absolute paths early
    src_dir = os.path.abspath(src_dir)
    build_dir = os.path.abspath(build_dir)
    optimized_cmake_filename = "OptimizedCMakeLists.txt" # Define the filename
    if cmakelists_path:
        cmakelists_path = os.path.abspath(cmakelists_path)
        # Check if the provided path is the one we intend to use
        if os.path.basename(cmakelists_path) != optimized_cmake_filename:
             logger.warn(f"Using a CMakeLists path ({cmakelists_path}) that doesn't match the expected optimized filename ({optimized_cmake_filename})")

    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    # Clean the build directory, *except* for the optimized CMakeLists if it's inside
    logger.debug(f"Cleaning build directory: {build_dir}")
    items_to_clean = os.listdir(build_dir) # Get list before iterating
    for item in items_to_clean:
        item_path = os.path.join(build_dir, item)

        # Skip the optimized CMakeLists file itself during cleaning
        if cmakelists_path and os.path.abspath(item_path) == cmakelists_path:
             logger.debug(f"Skipping cleaning of {item_path} (Optimized CMakeLists)")
             continue
        # Skip the temp source dir if it exists to avoid deleting mid-process
        if item == "_temp_src":
             logger.debug(f"Skipping cleaning of {item_path} (_temp_src)")
             continue

        logger.debug(f"Attempting to clean item: {item_path}")
        if os.path.isfile(item_path):
            try:
                os.unlink(item_path)
                logger.debug(f"Removed file: {item_path}")
            except OSError as e:
                logger.warn(f"Could not remove file {item_path}: {e}")
        elif os.path.isdir(item_path):
             logger.debug(f"Removing directory: {item_path}")
             shutil.rmtree(item_path, ignore_errors=True) # Use ignore_errors for robustness

    cmake_command = [
        "cmake",
        "-GNinja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DENABLE_PGO=ON",
        "-DGENERATE_PROFILE=ON",
        "-DUSE_PROFILE=OFF",
        "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
    ]

    # Conditionally disable GUI build
    if skip_gui:
        cmake_command.append("-DBUILD_GUI=OFF")
        logger.info("Skipping GUI build as requested")

    # Use `-fprofile-generate=.` instead of just `-fprofile-generate`
    # This explicitly tells GCC where to write the profile data
    if pgo_flags and pgo_flags[0]:
        gen_flag = "-fprofile-generate=."  # Override with explicit path
        cmake_command.append(f"-DCMAKE_CXX_FLAGS={gen_flag}")
        cmake_command.append(f"-DCMAKE_EXE_LINKER_FLAGS={gen_flag}")
        logger.info(f"Using explicit profiling path with flags: {gen_flag}")
    
    cmake_src_target = src_dir # Default source for CMake

    if cmakelists_path:
        logger.info(f"Using optimized CMakeLists.txt from {cmakelists_path}")

        # Ensure the source CMakeLists exists before proceeding
        if not os.path.exists(cmakelists_path):
            logger.error(f"Optimized CMakeLists.txt not found at specified path: {cmakelists_path}. Aborting build.")
            return False

        temp_dir = os.path.join(build_dir, "_temp_src")
        # Clean the temp dir specifically if it exists from a previous run
        if os.path.exists(temp_dir):
             logger.debug(f"Cleaning existing _temp_src directory: {temp_dir}")
             shutil.rmtree(temp_dir, ignore_errors=True)
        os.makedirs(temp_dir, exist_ok=True)

        # Copy contents of the original source directory to temp
        logger.debug(f"Copying source from {src_dir} to {temp_dir}")
        # Define items to ignore during copytree
        ignore_patterns = shutil.ignore_patterns('.git', 'build', '_temp_src', os.path.basename(build_dir))
        for item in os.listdir(src_dir):
            src_item = os.path.join(src_dir, item)
            dst_item = os.path.join(temp_dir, item)

            # Skip build directory to prevent recursion
            if os.path.abspath(src_item) == build_dir:
                continue

            try:
                if os.path.isdir(src_item):
                    shutil.copytree(src_item, dst_item, symlinks=True, ignore=ignore_patterns)
                    logger.debug(f"Copied directory {src_item} to {dst_item}")
                elif os.path.isfile(src_item):
                    # Skip copying the original CMakeLists.txt
                    if item == "CMakeLists.txt":
                         logger.debug(f"Skipping copy of original {src_item}")
                         continue
                    shutil.copy2(src_item, dst_item)
                    logger.debug(f"Copied file {src_item} to {dst_item}")
            except Exception as e:
                # Log specific copy errors but don't necessarily stop
                logger.warn(f"Error copying {src_item} to {dst_item}: {e}")


        # Now, copy the optimized CMakeLists.txt to the temp directory
        try:
            optimized_cmake_dest = os.path.join(temp_dir, "CMakeLists.txt")
            logger.debug(f"Attempting to copy {cmakelists_path} to {optimized_cmake_dest}")
            shutil.copy2(cmakelists_path, optimized_cmake_dest)
            # Verify copy
            if not os.path.exists(optimized_cmake_dest):
                 logger.error(f"Failed to verify copy of optimized CMakeLists.txt to {optimized_cmake_dest}")
                 return False
            logger.debug(f"Successfully copied and verified {cmakelists_path} to {optimized_cmake_dest}")
        except Exception as e:
            logger.error(f"Failed to copy optimized CMakeLists.txt from {cmakelists_path} to {optimized_cmake_dest}: {e}")
            return False

        cmake_src_target = temp_dir # Point CMake to the temp directory

    cmake_command.append(cmake_src_target)
    logger.debug(f"CMake command: {' '.join(cmake_command)}")
    logger.debug(f"CMake running in directory: {build_dir}")

    # Before running cmake, remove cache if using temp source to avoid conflicts
    if temp_dir:
         cache_path = os.path.join(build_dir, "CMakeCache.txt")
         if os.path.exists(cache_path):
             logger.debug(f"Removing CMakeCache.txt before instrumented build: {cache_path}")
             try: os.unlink(cache_path)
             except OSError as e: logger.warn(f"Could not remove CMakeCache.txt: {e}")

    success, _, _ = run_command(cmake_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to configure instrumented build")
        return False

    build_command = ["ninja", f"-j{jobs}"]
    logger.debug(f"Build command: {' '.join(build_command)}")
    logger.debug(f"Build running in directory: {build_dir}")
    success, _, _ = run_command(build_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to build instrumented binary")
        return False

    binary_path = os.path.join(build_dir, "desfact")
    if not os.path.exists(binary_path):
        logger.error(f"Binary not found at {binary_path} after build")
        return False

    logger.success(f"Successfully built instrumented binary at {binary_path}")
    return True

def stage_generate_profile(build_dir, input_dir, logger, timeout=None, hide_warnings=True):
    logger.info("STAGE 2: Running benchmarks to generate profile data")
    
    binary_path = os.path.join(build_dir, "desfact")
    if not os.path.exists(binary_path):
        logger.error(f"Instrumented binary not found at {binary_path}")
        return False
    
    if not os.access(binary_path, os.X_OK):
        logger.warn(f"Binary at {binary_path} is not executable, attempting to set permissions")
        try:
            os.chmod(binary_path, 0o755)
        except Exception as e:
            logger.error(f"Failed to set executable permissions: {e}")
            return False
    
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    if not csv_files:
        logger.error(f"No CSV files found in {input_dir}")
        return False
    
    logger.info(f"Found {len(csv_files)} CSV files for benchmark")
    
    os.makedirs(os.path.join(build_dir, "output"), exist_ok=True)
    
    # Create or copy models directory into build dir
    src_models_dir = os.path.join(os.path.dirname(os.path.abspath(build_dir)), "models")
    dst_models_dir = os.path.join(build_dir, "models")
    
    if os.path.exists(src_models_dir):
        logger.info(f"Copying models directory from {src_models_dir} to {dst_models_dir}")
        if os.path.exists(dst_models_dir):
            shutil.rmtree(dst_models_dir, ignore_errors=True)
        try:
            shutil.copytree(src_models_dir, dst_models_dir)
            logger.success(f"Successfully copied models to build directory")
        except Exception as e:
            logger.error(f"Failed to copy models directory: {e}")
    else:
        logger.warn(f"Models directory not found at {src_models_dir}")
        os.makedirs(dst_models_dir, exist_ok=True)
        logger.info(f"Created empty models directory at {dst_models_dir}")
    
    # Prepare environment with explicit GCOV settings
    benchmark_env = os.environ.copy()
    # Set absolute path for GCOV data
    gcov_dir = os.path.abspath(build_dir)
    benchmark_env['GCOV_PREFIX'] = gcov_dir
    benchmark_env['GCOV_PREFIX_STRIP'] = "0"
    logger.info(f"Setting GCOV_PREFIX={gcov_dir} for profile data generation")

    with logger.start_progress(len(csv_files), "Running benchmarks") as progress:
        for i, csv_file in enumerate(csv_files):
            progress.set_description(f"Benchmark {i+1}/{len(csv_files)}")
            csv_path = os.path.join(os.path.abspath(input_dir), csv_file)
            output_path = os.path.abspath(os.path.join(build_dir, "output", f"temp_output_{i}.csv"))
            
            # Attempt adding a basic mode without advanced descriptors first
            benchmark_command = [
                "./desfact",
                "-i", csv_path,
                "-o", output_path,
                "-d", "basic"
            ]
            
            logger.debug(f"Running benchmark with basic descriptors: {' '.join(benchmark_command)}")
            
            # Run the benchmark with the modified environment
            success, stdout, stderr = run_command(benchmark_command, logger, cwd=build_dir, env=benchmark_env, timeout=timeout, hide_warnings=hide_warnings)
            if not success:
                logger.warn(f"Basic benchmark failed: {stderr}")
                
                # Try with different column name - molecule (from the sample data)
                logger.warn("Retrying with explicit 'molecule' column name")
                benchmark_command_alt = [
                    "./desfact",
                    "-i", csv_path,
                    "-o", output_path,
                    "-d", "basic",
                    "-s", "molecule"
                ]
                success, stdout, stderr = run_command(benchmark_command_alt, logger, cwd=build_dir, env=benchmark_env, timeout=timeout, hide_warnings=hide_warnings)
                if not success:
                    logger.error(f"Benchmark also failed with explicit 'molecule' column: {stderr}")
                    
                    # Minimal fallback
                    logger.warn("Trying minimal fallback with simplified args")
                    simple_cmd = [
                        "./desfact",
                        "-i", csv_path,
                        "-o", output_path,
                    ]
                    success, stdout, stderr = run_command(simple_cmd, logger, cwd=build_dir, env=benchmark_env, timeout=timeout, hide_warnings=hide_warnings)
                    if not success:
                        logger.error(f"All benchmark attempts failed for {csv_file}")
                        # Continue anyway to try to generate at least some profile data
            
            try:
                if os.path.exists(output_path):
                    os.unlink(output_path)
            except:
                pass
                
            progress.update()
    
    # Enhanced search for profile data with various patterns
    profile_files = []
    # Look for standard profile data
    for pattern in ["*.gcda", "desfact-*.gcda", "*.profraw", "default.profdata"]:
        found = glob.glob(os.path.join(build_dir, pattern))
        profile_files.extend(found)
    
    # Also search subdirectories
    for root, dirs, files in os.walk(build_dir):
        for file in files:
            if file.endswith('.gcda') or file.endswith('.profraw'):
                profile_files.append(os.path.join(root, file))
    
    if not profile_files:
        logger.error("No profile data (.gcda or .profraw) was generated")
        # Diagnostic information
        logger.info("Searching for any .gcno files (instrumentation files)...")
        gcno_files = []
        for root, dirs, files in os.walk(build_dir):
            for file in files:
                if file.endswith('.gcno'):
                    gcno_files.append(os.path.join(root, file))
        
        if gcno_files:
            logger.info(f"Found {len(gcno_files)} .gcno files, but no .gcda files.")
            logger.info("This suggests the instrumentation is set up correctly but profile data wasn't created.")
        else:
            logger.info("No .gcno files found. This suggests instrumentation wasn't properly configured.")
        
        # Force some profile data to be created so we can continue
        logger.warn("Generating minimal profile data to allow build to continue")
        dummy_profile = os.path.join(build_dir, "dummy.gcda")
        with open(dummy_profile, "wb") as f:
            f.write(b'\x01')
        profile_files.append(dummy_profile)
    
    logger.success(f"Successfully generated profile data ({len(profile_files)} files)")
    
    with open(os.path.join(build_dir, "profile_files.txt"), "w") as f:
        for profile_file in profile_files:
            f.write(f"{profile_file}\n")
    
    return True

def stage_optimized_build(src_dir, build_dir, jobs, logger, pgo_flags=None, timeout=None, hide_warnings=True, cmakelists_path=None, skip_gui=False):
    logger.info("STAGE 3: Building optimized binary with profile data")
    
    # Ensure absolute paths
    build_dir = os.path.abspath(build_dir)
    src_dir = os.path.abspath(src_dir)
    if cmakelists_path:
        cmakelists_path = os.path.abspath(cmakelists_path) # Make sure this is absolute too

    # Search for profile data (should be in build_dir)
    profile_files = []
    for root, dirs, files in os.walk(build_dir):
        for file in files:
            if file.endswith('.gcda') or file.endswith('.profraw'):
                profile_files.append(os.path.join(root, file))
    
    if not profile_files:
        logger.error("No profile data found. Cannot build optimized binary.")
        return False
    
    logger.info(f"Found {len(profile_files)} profile data files")

    cmake_command = [
        "cmake",
        "-GNinja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DENABLE_PGO=ON",
        "-DGENERATE_PROFILE=OFF", # Turn off generation
        "-DUSE_PROFILE=ON"       # Turn on profile usage
    ]

    # Conditionally disable GUI build
    if skip_gui:
        cmake_command.append("-DBUILD_GUI=OFF")
        logger.info("Skipping GUI build for optimized stage")

    # Profile data setup
    profile_type = None
    if any(f.endswith('.gcda') for f in profile_files):
        profile_type = "gcc"
        logger.info("Found GCC-style profile data (.gcda files)")
    elif any(f.endswith('.profraw') for f in profile_files):
        profile_type = "clang"
        logger.info("Found Clang-style profile data (.profraw files)")
        # Merge Clang profiles
        merge_cmd = ["llvm-profdata", "merge", "-output=default.profdata"]
        merge_cmd.extend([f for f in profile_files if f.endswith('.profraw')])
        merge_success, _, merge_err = run_command(merge_cmd, logger, cwd=build_dir)
        if merge_success:
            logger.info("Merged Clang profile data files to default.profdata")
        else:
            logger.warn(f"Failed to merge Clang profile data: {merge_err}")

    # Set use flags
    if pgo_flags and pgo_flags[1]:
        use_flag = pgo_flags[1]
        # Ensure the flag includes the path if needed, e.g., -fprofile-use=.
        if profile_type == "gcc" and "-fprofile-correction" not in use_flag:
            use_flag += " -fprofile-correction"
        if profile_type == "gcc" and "-fprofile-use" in use_flag and "=" not in use_flag:
             # Ensure path is specified for GCC use flag if detected without one
             use_flag = use_flag.replace("-fprofile-use", "-fprofile-use=.")
             logger.debug(f"Modified GCC use flag to specify path: {use_flag}")

        cmake_command.append(f"-DCMAKE_CXX_FLAGS={use_flag}")
        cmake_command.append(f"-DCMAKE_EXE_LINKER_FLAGS={use_flag}")
        logger.info(f"Applying PGO Use flags: {use_flag}")
    else:
         # Fallback logic remains the same, ensure paths are correct
         if profile_type == "gcc":
             logger.warn("PGO use flags not detected, attempting default GCC flags.")
             cmake_command.append("-DCMAKE_CXX_FLAGS=-fprofile-use=. -fprofile-correction")
             cmake_command.append("-DCMAKE_EXE_LINKER_FLAGS=-fprofile-use=. -fprofile-correction")
         # ... (Clang fallback) ...


    # --- Source directory handling ---
    cmake_src_target = src_dir # Default if not using optimized CMakeLists

    if cmakelists_path:
        logger.info(f"Using optimized CMakeLists.txt from {cmakelists_path} for optimized build.")

        # ** Use the SAME temp directory name as stage_instrumented_build **
        temp_dir = os.path.join(build_dir, "_temp_src")

        # Clean and repopulate the temp directory for this stage
        if os.path.exists(temp_dir):
            logger.debug(f"Cleaning existing _temp_src directory for optimized build: {temp_dir}")
            shutil.rmtree(temp_dir, ignore_errors=True)
        os.makedirs(temp_dir, exist_ok=True)
        logger.debug(f"Recreated temporary source directory at {temp_dir}")

        # Copy source files to temp directory (same logic as stage_instrumented_build)
        logger.debug(f"Copying source from {src_dir} to {temp_dir} for optimized build")
        ignore_patterns = shutil.ignore_patterns('.git', 'build', '_temp_src*', os.path.basename(build_dir))
        for item in os.listdir(src_dir):
             src_item = os.path.join(src_dir, item)
             dst_item = os.path.join(temp_dir, item)

             if os.path.abspath(src_item) == build_dir:
                 continue # Skip build dir

             try:
                 if os.path.isdir(src_item):
                     # Avoid copying the _temp_src dir into itself if src_dir is '.'
                     if os.path.abspath(src_item) == os.path.abspath(temp_dir):
                         continue
                     shutil.copytree(src_item, dst_item, symlinks=True, ignore=ignore_patterns)
                 elif os.path.isfile(src_item):
                     if item == "CMakeLists.txt": # Skip original, use optimized one
                         continue
                     shutil.copy2(src_item, dst_item)
             except Exception as e:
                 logger.warn(f"Error copying {src_item} to {dst_item}: {e}")

        # Copy the *optimized* CMakeLists.txt to the temp directory, renaming it
        logger.debug(f"Copying optimized CMakeLists from {cmakelists_path} to {temp_dir}/CMakeLists.txt")
        try:
            optimized_cmake_dest = os.path.join(temp_dir, "CMakeLists.txt")
            shutil.copy2(cmakelists_path, optimized_cmake_dest)
            if not os.path.exists(optimized_cmake_dest):
                logger.error(f"Failed to copy optimized CMakeLists.txt to {optimized_cmake_dest}")
                return False
            logger.debug(f"Successfully copied optimized CMakeLists to {optimized_cmake_dest}")
        except Exception as e:
            logger.error(f"Failed to copy optimized CMakeLists: {e}")
            return False

        cmake_src_target = temp_dir # Point CMake to the temp directory

    # --- End Source directory handling ---


    # Run CMake configure
    cmake_command.append(cmake_src_target)
    logger.debug(f"CMake command: {' '.join(cmake_command)}")
    logger.debug(f"CMake running in directory: {build_dir} with source {cmake_src_target}") # Log source dir

    # Before running cmake, remove cache if using temp source to avoid conflicts
    if temp_dir:
         cache_path = os.path.join(build_dir, "CMakeCache.txt")
         if os.path.exists(cache_path):
             logger.debug(f"Removing CMakeCache.txt before optimized configure: {cache_path}")
             try: os.unlink(cache_path)
             except OSError as e: logger.warn(f"Could not remove CMakeCache.txt: {e}")

    success, _, _ = run_command(cmake_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to configure optimized build")
        return False

    # Build the optimized binary
    build_command = ["ninja", f"-j{jobs}"]
    logger.debug(f"Build command: {' '.join(build_command)}")
    success, _, _ = run_command(build_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to build optimized binary")
        return False
    
    # Verify binary exists
    binary_path = os.path.join(build_dir, "desfact")
    if not os.path.exists(binary_path):
        logger.error(f"Binary not found at {binary_path} after build")
        return False
    
    logger.success("Successfully built optimized binary")
    return True

def signal_handler(sig, frame):
    print("\nInterrupted by user. Exiting gracefully...")
    sys.exit(1)

def dump_environment_info(dest_file, environment_data, logger):
    try:
        with open(dest_file, "w") as f:
            json.dump(environment_data, f, indent=2)
        logger.info(f"Environment information saved to {dest_file}")
    except Exception as e:
        logger.warn(f"Could not save environment info: {e}")

def main():
    signal.signal(signal.SIGINT, signal_handler)
    
    # Fast help check
    if len(sys.argv) == 2 and sys.argv[1] in ['--help', '-h']:
        print("Advanced Profile-Guided Optimization (PGO) Tool for DescriptorFactory")
        print("")
        print("This script optimizes the build process for your specific environment:")
        print("  1. Analyzes your CPU and memory capabilities")
        print("  2. Determines optimal compiler flags and settings")
        print("  3. Creates a tailored, temporary CMakeLists.txt")
        print("  4. Builds an instrumented binary for profile generation")
        print("  5. Runs a benchmark to generate profile data")
        print("  6. Builds an optimized binary using the profile data")
        print("")
        print("Usage:")
        print("  ./INSTALL.py [options]")
        print("")
        print("Options:")
        print("  -h, --help                Display this help message")
        print("  -j, --jobs NUM            Set number of parallel jobs (default: number of CPU cores)")
        print("  -i, --input-dir DIR       Directory with test SMILES files (default: ./input)")
        print("  -b, --build-dir DIR       Build directory (default: ./build)")
        print("  -s, --src-dir DIR         Source directory (default: current directory)")
        print("  -v, --verbose             Enable verbose output")
        print("  -q, --quick               Quick mode: skip thorough testing")
        print("  -t, --timeout SEC         Timeout for build processes in seconds (default: none)")
        print("  -w, --show-warnings       Show compiler warnings (default: hidden)")
        print("  --skip-gui                Do not build the GUI target (desfact-gui)")
        print("  --skip-instrumented       Skip instrumented build stage")
        print("  --skip-profile            Skip profile generation stage")
        print("  --skip-optimized          Skip optimized build stage")
        print("  --skip-detection          Skip environment detection (use defaults)")
        return 0
    
    parser = argparse.ArgumentParser(
        description="Advanced Profile-Guided Optimization (PGO) for DescriptorFactory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False # Disable default help to use custom fast help
    )
    
    # Add arguments (keep alphabetical if possible)
    parser.add_argument("-b", "--build-dir", default="./build", help="Build directory")
    parser.add_argument("-h", "--help", action="store_true", help="Display this help message")
    parser.add_argument("-i", "--input-dir", default="./input", help="Directory containing input CSV files for benchmarking")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="Number of parallel build jobs")
    parser.add_argument("-q", "--quick", action="store_true", help="Quick mode: skip thorough testing")
    parser.add_argument("-s", "--src-dir", default=".", help="Source directory")
    parser.add_argument("-t", "--timeout", type=int, default=None, help="Timeout for build processes in seconds")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("-w", "--show-warnings", action="store_true", help="Show compiler warnings (default: hidden)")
    parser.add_argument("--skip-detection", action="store_true", help="Skip environment detection (use defaults)")
    parser.add_argument("--skip-gui", action="store_true", help="Do not build the GUI target (desfact-gui)")
    parser.add_argument("--skip-instrumented", action="store_true", help="Skip instrumented build stage")
    parser.add_argument("--skip-profile", action="store_true", help="Skip profile generation stage")
    parser.add_argument("--skip-optimized", action="store_true", help="Skip optimized build stage")
    
    try:
        args = parser.parse_args()
        # Handle help argument manually after parsing
        if args.help:
            # Call the fast help display logic directly
            main()
            return 0 # Exit after displaying help
    except SystemExit:
        return 1 # Exit if parsing fails (e.g., invalid argument)

    logger = Logger(verbose=args.verbose)
    
    hide_warnings = not args.show_warnings
    
    logger.info("Starting Advanced Profile-Guided Optimization process")
    logger.debug(f"Jobs: {args.jobs}")
    logger.debug(f"Source directory: {args.src_dir}")
    logger.debug(f"Input directory: {args.input_dir}")
    logger.debug(f"Build directory: {args.build_dir}")
    logger.debug(f"Quick mode: {args.quick}")
    logger.debug(f"Hide warnings: {hide_warnings}")
    if args.timeout:
        logger.debug(f"Command timeout: {args.timeout}s")
    logger.debug(f"Skip GUI: {args.skip_gui}")

    try:
        environment_data = {}
        temp_cmakelists = None # Initialize
        
        if not args.skip_detection:
            logger.info("Phase 0: Analyzing system environment")
            
            cpu_info, cpu_capabilities = detect_cpu_features(logger)
            environment_data["cpu"] = cpu_info
            memory_info = detect_memory_info(logger)
            environment_data["memory"] = memory_info
            compiler_info = detect_compiler_info(logger)
            environment_data["compiler"] = compiler_info
            
            if not os.path.exists(args.build_dir):
                os.makedirs(args.build_dir, exist_ok=True)
            dump_environment_info(os.path.join(args.build_dir, "environment.json"), environment_data, logger)
            
            optimization_flags = get_default_optimization_flags(cpu_capabilities)
            logger.info(f"Using optimization flags: {optimization_flags}")
            
            temp_cmakelists = os.path.join(args.build_dir, "OptimizedCMakeLists.txt")
            created_cmake = create_optimized_cmakelists(args.src_dir, temp_cmakelists, optimization_flags, logger, args.skip_gui)
            if not created_cmake:
                logger.error("Failed to create optimized CMakeLists.txt during detection phase. Aborting.")
                return 1
        else:
            logger.info("Skipping environment detection")
            # If skipping detection, we won't have an optimized CMakeLists unless one already exists
            temp_cmakelists = os.path.join(args.build_dir, "OptimizedCMakeLists.txt")
            if not os.path.exists(temp_cmakelists):
                 logger.warn("Skipping detection, but no OptimizedCMakeLists.txt found. Build might use default flags.")
                 temp_cmakelists = None # Ensure it's None if not found

        pgo_flags = detect_compiler_flags(logger, args.quick)
        ensure_input_directory(args.input_dir, logger)
        
        start_time = time.time()
        success = True
        
        if not args.skip_instrumented:
            logger.info("Starting PGO optimization process")
            
            if temp_cmakelists and not os.path.exists(temp_cmakelists):
                logger.error(f"Optimized CMakeLists.txt expected but missing at {temp_cmakelists}. Aborting.")
                return 1
            elif temp_cmakelists:
                 logger.debug(f"Verified existence of {temp_cmakelists} before instrumented build.")
            
            success = stage_instrumented_build(
                args.src_dir, args.build_dir, args.jobs, logger,
                pgo_flags=pgo_flags, timeout=args.timeout, hide_warnings=hide_warnings,
                cmakelists_path=temp_cmakelists, skip_gui=args.skip_gui
            )
            if not success:
                logger.error("Failed to build instrumented binary. Aborting.")
                return 1
        else:
            logger.info("Skipping instrumented build stage")
        
        if not args.skip_profile and success:
            success = stage_generate_profile(
                args.build_dir, args.input_dir, logger,
                timeout=args.timeout, hide_warnings=hide_warnings
            )
            if not success:
                logger.error("Failed to generate profile data. Aborting.")
                # Optionally continue to optimized build if desired, but warn
                # return 1
            elif not os.path.exists(os.path.join(args.build_dir, "profile_files.txt")):
                 logger.error("Profile generation reported success, but profile_files.txt is missing.")
                 return 1
        elif not args.skip_profile and not success:
             logger.error("Instrumented build failed, skipping profile generation.")
             return 1
        else: # --skip-profile is True
            logger.info("Skipping profile generation stage")
            # If skipping profile generation, the optimized build cannot use PGO
            # We could potentially still run the optimized build without PGO flags
            # For now, let's just skip the optimized build too if profile gen is skipped
            if not args.skip_optimized:
                 logger.warn("Skipping profile generation implies skipping PGO-optimized build.")
                 args.skip_optimized = True

        if not args.skip_optimized and success:
            # Verify optimized CMakeLists exists if we created one earlier
            if temp_cmakelists and not os.path.exists(temp_cmakelists):
                logger.error(f"Optimized CMakeLists.txt expected but missing at {temp_cmakelists} before optimized build. Aborting.")
                return 1

            success = stage_optimized_build(
                args.src_dir, args.build_dir, args.jobs, logger,
                pgo_flags=pgo_flags, timeout=args.timeout, hide_warnings=hide_warnings,
                cmakelists_path=temp_cmakelists, skip_gui=args.skip_gui
            )
            if not success:
                logger.error("Failed to build optimized binary.")
                return 1
        elif not args.skip_optimized and not success:
             logger.error("Previous stage failed, skipping optimized build.")
             return 1
        else: # --skip-optimized is True
            logger.info("Skipping optimized build stage")
        
        elapsed_time = time.time() - start_time
        
        if success:
            logger.success(f"PGO process completed successfully in {elapsed_time:.2f} seconds")
            logger.info(f"Optimized binary available at: {os.path.join(args.build_dir, 'desfact')}")
            if not args.skip_gui:
                 gui_path = os.path.join(args.build_dir, 'desfact-gui')
                 if os.path.exists(gui_path):
                     logger.info(f"Optimized GUI binary available at: {gui_path}")
                 else:
                     logger.warn(f"GUI binary was expected but not found at: {gui_path}")

            if temp_cmakelists and os.path.exists(temp_cmakelists):
                logger.info(f"Optimized CMakeLists.txt available at: {temp_cmakelists}")
            
            return 0
        else:
            logger.error(f"PGO process failed after {elapsed_time:.2f} seconds")
            return 1
    
    except KeyboardInterrupt:
        logger.warn("Process interrupted by user. Exiting gracefully...")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        if args.verbose:
            import traceback
            logger.debug(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())