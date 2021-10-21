#include "tests/testing.hpp"
#include "tests/test_penning_trap.hpp"
#include "tests/test_runge_kutta_forward_euler.hpp"

void run_testing() {
    test_penning_trap();
    test_runge_kutta_forward_euler();
}
