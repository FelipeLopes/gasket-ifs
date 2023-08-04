cc_library(
    name = "gasket",
    srcs = glob(["src/gasket/*.cpp"]),
    hdrs = glob(["src/gasket/*.hpp"]),
    deps = [
        "@gmp//:gmp",
        "@boost//:boost",
    ],
)

cc_binary(
    name = "gasket-ifs",
    srcs = ["src/main.cpp"],
    deps = [
        ":gasket",
    ],
)