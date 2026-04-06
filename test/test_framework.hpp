#ifndef TEST_FRAMEWORK_HPP
#define TEST_FRAMEWORK_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>

namespace test {

struct TestResult {
    std::string name;
    bool passed;
    std::string message;
};

class TestRunner {
public:
    static TestRunner& instance() {
        static TestRunner runner;
        return runner;
    }

    void add_test(const std::string& name, std::function<bool()> test_func) {
        tests_.push_back({name, test_func});
    }

    int run() {
        int passed = 0;
        int failed = 0;

        std::cout << "\n========================================\n";
        std::cout << "Running " << tests_.size() << " tests...\n";
        std::cout << "========================================\n\n";

        for (const auto& [name, func] : tests_) {
            std::cout << "[TEST] " << name << " ... ";
            try {
                if (func()) {
                    std::cout << "✓ PASS\n";
                    ++passed;
                } else {
                    std::cout << "✗ FAIL\n";
                    ++failed;
                }
            } catch (const std::exception& e) {
                std::cout << "✗ EXCEPTION: " << e.what() << "\n";
                ++failed;
            }
        }

        std::cout << "\n========================================\n";
        std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
        std::cout << "========================================\n";

        return failed;
    }

private:
    std::vector<std::pair<std::string, std::function<bool()>>> tests_;
};

#define TEST(name) \
    static bool test_##name(); \
    static struct test_##name##_register { \
        test_##name##_register() { \
            test::TestRunner::instance().add_test(#name, test_##name); \
        } \
    } test_##name##_instance; \
    static bool test_##name()

#define ASSERT_TRUE(expr) \
    if (!(expr)) { \
        std::cerr << "\n  Assertion failed: " #expr " (expected true, got false)\n"; \
        return false; \
    }

#define ASSERT_FALSE(expr) \
    if (expr) { \
        std::cerr << "\n  Assertion failed: " #expr " (expected false, got true)\n"; \
        return false; \
    }

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::cerr << "\n  Assertion failed: " #a " == " #b "\n"; \
        std::cerr << "  Left: " << (a) << "\n"; \
        std::cerr << "  Right: " << (b) << "\n"; \
        return false; \
    }

#define ASSERT_NE(a, b) \
    if ((a) == (b)) { \
        std::cerr << "\n  Assertion failed: " #a " != " #b "\n"; \
        return false; \
    }

#define ASSERT_NEAR(a, b, epsilon) \
    if (std::abs((a) - (b)) > (epsilon)) { \
        std::cerr << "\n  Assertion failed: |" << (a) << " - " << (b) << "| <= " << (epsilon) << "\n"; \
        return false; \
    }

#define ASSERT_THROW(expr, exc_type) \
    { \
        bool caught = false; \
        try { \
            expr; \
        } catch (const exc_type&) { \
            caught = true; \
        } catch (...) { \
            std::cerr << "\n  Assertion failed: " #expr " threw wrong exception type\n"; \
            return false; \
        } \
        if (!caught) { \
            std::cerr << "\n  Assertion failed: " #expr " did not throw " #exc_type "\n"; \
            return false; \
        } \
    }

#define ASSERT_NO_THROW(expr) \
    try { \
        expr; \
    } catch (const std::exception& e) { \
        std::cerr << "\n  Assertion failed: " #expr " threw exception: " << e.what() << "\n"; \
        return false; \
    }

#define RUN_TESTS() \
    return test::TestRunner::instance().run()

} // namespace test

#endif // TEST_FRAMEWORK_HPP
