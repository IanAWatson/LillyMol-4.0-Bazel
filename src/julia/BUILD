cc_library (
    name = "lillymol_julia",
    srcs = [
        "lillymol.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "@libcxxwrap//:libcxxwrap",
        "@julia//:julia",
    ],
)

cc_shared_library(
    name = "lillymol_shared",
    deps = [
        "lillymol_julia",
    ],
    dynamic_deps = [
    ],
    shared_lib_name = "lillymol_julia.so",
)
