import os
from setuptools import setup, Extension
from distutils.sysconfig import customize_compiler
import Cython.Build
from Cython.Build import cythonize

def CppExtension(name, sources):
    
    include_dirs=[
            "src/",
            "external/wfa2/", 
            "external/wfa2/include/", 
            "external/wfa2/bindings/cpp/",
            "external/wrappers/"
        ]
    library_dirs=["external/wfa2/lib/"]
    libraries=["wfacpp"]
    
    return Extension(
        name,
        sources=sources,
        language="c++",
        extra_compile_args=["-std=c++20", "-Werror=return-type", "-Werror=narrowing"],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=libraries,
        undef_macros=["NDEBUG"],
    )


extensions = [
    CppExtension(
        "giggles.core",
        sources=[
            "giggles/core.pyx",
            "src/binomial.cpp",
            "src/bipartitioniterator.cpp",
            "src/column.cpp",
            "src/columniterator.cpp",
            "src/emissionprobabilitycomputer.cpp",
            "src/entry.cpp",
            "src/genotype.cpp",
            "src/genotypehmm.cpp",
            "src/genotypelikelihoods.cpp",
            "src/genotypingalgorithm.cpp",
            "src/graycodes.cpp",
            "src/haplotypemapper.cpp",
            "src/indexset.cpp",
            "src/read.cpp",
            "src/readset.cpp",
            "src/transitionprobabilitycomputer.cpp",
            "src/phasing/phasingcolumncostcomputer.cpp",
            "src/phasing/phasingcolumnindexingiterator.cpp",
            "src/phasing/phasingcolumnindexingscheme.cpp",
            "src/phasing/phasingcolumniterator.cpp",
            "src/phasing/phasingdptable.cpp",
        ]
    ),
    CppExtension("giggles.readselect", sources=["giggles/readselect.pyx"]),
    CppExtension("giggles.cpptests", sources=["giggles/cpptests.pyx"]),
    CppExtension(
        "giggles.align", 
        sources=[
            "giggles/align.pyx",
            "external/wrappers/wfawrapper.cpp"
        ]
    ),
    CppExtension("giggles.priorityqueue", sources=["giggles/priorityqueue.pyx"]),
    CppExtension("giggles._variants", sources=["giggles/_variants.pyx"]),
]


class BuildExt(Cython.Build.build_ext):
    def build_extensions(self):
        # Remove the warning about “-Wstrict-prototypes” not being valid for C++,
        # see http://stackoverflow.com/a/36293331/715090
        customize_compiler(self.compiler)
        if self.compiler.compiler_so[0].endswith("clang"):
            # Clang needs this option in order to find the unordered_set header
            print("detected clang, using option -stdlib=libc++")
            self.compiler.compiler_so.append("-stdlib=libc++")
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        super().build_extensions()


# Avoid compilation if we are being installed within Read The Docs
if os.environ.get("READTHEDOCS") == "True":
    cmdclass = {}
    ext_modules = []
    install_requires = []
else:
    cmdclass = {"build_ext": BuildExt}
    ext_modules = extensions
    install_requires = [
        "pysam>=0.18.0",
        "pyfaidx>=0.5.5.2",
        "biopython>=1.73",  # pyfaidx needs this for reading bgzipped FASTA files
        "xopen>=1.2.0"
    ]

setup(
    use_scm_version={"write_to": "giggles/_version.py"},
    cmdclass=cmdclass,
    ext_modules=cythonize(ext_modules),
    install_requires=install_requires,
)
