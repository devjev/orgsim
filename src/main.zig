const std = @import("std");
const testing = std.testing;
const trait = std.meta.trait;

pub const GammaFunctionError = error{
    ArgumentMustBeGreaterThanZero,
};

// Natural logarithm of the Gamma function.
// Requires a floating type parameter. Should be at least f32 or bigger.
pub fn log_gamma(comptime float_t: type, x: float_t) GammaFunctionError!float_t {
    comptime if (!trait.isFloat(float_t) or @sizeOf(float_t) < @sizeOf(f32)) {
        @compileError("log_gamma requires at least a 32-bit float type parameter.");
    };

    if (x <= 0.0) return GammaFunctionError.ArgumentMustBeGreaterThanZero;

    const coefs = [14]float_t{ 57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, 0.339946499848118887E-4, 0.465236289270485756E-4, -0.983744753048795646E-4, 0.158088703224912494E-3, -0.210264441724104883E-3, 0.217439618115212643E-3, -0.164318106536763890E-3, 0.844182239838527433E-4, -0.261908384015814087E-4, 0.368991826595316234E-5 };

    const stp: float_t = 2.5066282746310005;
    var ser: float_t = 0.999999999999997092;
    var yy: float_t = x;
    var xx: float_t = x;

    var tmp: float_t = xx + 5.2421875;

    tmp = (xx + 0.5) * @log(tmp) - tmp;

    inline for (coefs) |coef| {
        yy += 1.0;
        ser += coef / yy;
    }

    return tmp + @log(stp * ser / xx);
}

fn gser(comptime float_t: type, a: float_t, b: float_t) GammaFunctionError!float_t {
    _ = b;
    _ = a;

}

test "log_gamma(x + 1) = factorial(x) " {
    const float_t = f32;
    const test_tuple_t = struct {
        x: float_t,
        y: float_t,
    };

    const cases = [_]test_tuple_t{
        test_tuple_t{ .x = 2.0, .y = 2.0 },
        test_tuple_t{ .x = 3.0, .y = 6.0 },
        test_tuple_t{ .x = 4.0, .y = 24.0 },
        test_tuple_t{ .x = 5.0, .y = 120.0 },
        test_tuple_t{ .x = 6.0, .y = 720.0 },
        test_tuple_t{ .x = 7.0, .y = 5_040.0 },
        test_tuple_t{ .x = 8.0, .y = 40_320.0 },
        test_tuple_t{ .x = 9.0, .y = 362_880.0 },
        test_tuple_t{ .x = 10.0, .y = 3_628_800.0 },
    };

    const epsilon: float_t = 1E-5;

    inline for (cases) |case| {
        const expected_result: float_t = case.y;
        const log_got: float_t = try log_gamma(float_t, case.x + 1.0);
        const got = @exp(log_got);
        try testing.expectApproxEqRel(expected_result, got, epsilon);
    }
}
