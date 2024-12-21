# Followed https://dnmtechs.com/using-cmake-in-setup-py-extending-setuptools-for-python-3/
# with changes for customization
import sys
import os
import subprocess
import re
import shutil

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.module_name = name
        self.sourcedir = os.path.abspath(sourcedir)
    
class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cfg = 'Debug' if self.debug else 'Release'
        if is_editable_install:
            editable_opt = "ON"
        else:
            editable_opt = "OFF"
            extdir = os.path.join(extdir, ext.name)
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            f"-DSETUP_PY_BUILD=ON",
            f"-DCPP_BPOTF_VER={compose_version()}",
            f"-DEDITABLE_INST={editable_opt}"
        ]

        build_args = ['--config', cfg]

        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        self.add_platform_specific_args(cmake_args, build_args, cmake_generator, cfg, extdir)
        
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            if hasattr(self, "parallel") and self.parallel:
                build_args += [f"-j{self.parallel}"]

        # Configure and build the project
        build_temp = self.build_temp
        os.makedirs(build_temp, exist_ok=True)
        # print(cmake_args)
        subprocess.check_call(["cmake", os.path.abspath(".")] + cmake_args, cwd=build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)

    def add_platform_specific_args(self, cmake_args, build_args, cmake_generator, cfg, extdir):
        """
        Add platform-specific arguments to cmake_args and build_args.
        """
        if sys.platform.startswith("win"):
            self._add_windows_args(cmake_args, build_args, cmake_generator, cfg, extdir)
        elif sys.platform.startswith("darwin"):
            self._add_macos_args(cmake_args, build_args)
        elif sys.platform.startswith("linux"):
            self._add_linux_args(cmake_args, build_args)
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}")

    def _add_windows_args(self, cmake_args, build_args, cmake_generator, cfg, extdir):
        pass

    def _add_macos_args(self, cmake_args, build_args):
        # Handle macOS-specific flags
        archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
        if archs:
            cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]
        else:
            import platform
            arch = platform.machine()
            if arch:
                cmake_args += [f"-DCMAKE_OSX_ARCHITECTURES={arch}"]

    def _add_linux_args(self, cmake_args, build_args):
        # Handle Linux-specific flags
        if shutil.which("gcc-10") and shutil.which("g++-10"):
            os.environ["CC"] = "gcc-10"
            os.environ["CXX"] = "g++-10"

    def run(self):
        # Run the build process
        super().run()

        if not is_editable_install:
            module_name = self.extensions[0].module_name
            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(module_name)))
            extdir = os.path.join(extdir, module_name)
            os.makedirs(extdir, exist_ok=True)

            # Generate __init__.py in the output folder
            init_file = os.path.join(extdir, "__init__.py")
            with open(init_file, "w") as f:
                f.write("# This file is generated during installation\n")
                f.write("from .BPOTF import *\n")
                f.write("from .BPOTF import __version__\n")
                f.write("__all__ = ['__version__']\n")
    
class CustomInstallCommand(install):
    def run(self):
        super().run()

        self.generate_stubs()

    def generate_stubs(self):
        print("INFO -- Generating stubs!")
        build_ext_cmd = self.get_finalized_command("build_ext")
        module_name = build_ext_cmd.extensions[0].module_name
        build_temp = build_ext_cmd.build_temp

        env = os.environ.copy()
        env['PYTHONPATH'] = self.build_lib + os.pathsep + env.get('PYTHONPATH', '')
        
        stubgen_cmd = ["pybind11-stubgen", 
                       "--enum-class-locations", "NoiseType:BPOTF.BPOTF.OBPOTF.NoiseType", 
                       "-o"]
        
        # extdir = os.path.abspath(os.path.dirname(build_ext_cmd.get_ext_fullpath(module_name)))
        stubs_out_dir = os.path.join(build_temp, "stubs")

        subprocess.check_call(stubgen_cmd + [stubs_out_dir, module_name], env=env)

        bpotf_stub_filepath = os.path.join(stubs_out_dir, module_name, "BPOTF.pyi")
        dest_path = os.path.join(self.install_lib, module_name, "BPOTF.pyi")
        shutil.copyfile(bpotf_stub_filepath, dest_path)
        print(f"INFO -- Stub files copied to {dest_path}")
        

def compose_version() -> str:
    version = "dev0"
    with open("version", "r") as ver_file:
        version_components = [line.split()[1] for line in ver_file if len(line.split()) >= 2]
        version = ".".join(version_components)
    return version

os.makedirs("./build", exist_ok=True)

with open("README.md", 'r', encoding='utf-8') as file:
    long_des = file.read()

is_editable_install = "editable_wheel" in sys.argv

setup(
    name="BPOTF",
    version=compose_version(),
    author="Imanol Etxezarreta",
    author_email="ietxezarretam@gmail.com",
    url="https://github.com/Ademartio/BPOTF",
    description="Implementation of Belief Propagation Ordered Tanner Forest decoding method.",
    long_description=long_des,
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("BPOTF", "src")],
    cmdclass={
        "build_ext": CMakeBuild,
        "install": CustomInstallCommand},
    packages=["BPOTF"],
    zip_safe=False,
    python_requires=">=3.8",
)
