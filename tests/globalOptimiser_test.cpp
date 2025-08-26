#include "../src/GlobalOptimiser/globalOptimiser.hpp"
#include <gtest/gtest.h>
#include <Eigen/Dense>

// Mock classes for testing
class MockDesignVariable {
public:
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    
    MockDesignVariable() : vec(Vector::Zero(3)) {}
    explicit MockDesignVariable(const Vector& v) : vec(v) {}
    
    Vector& getVector() { return vec; }
    const Vector& getVector() const { return vec; }
private:
    Vector vec;
};

class MockResponse {
public:
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    
    MockResponse() : objectives(Vector::Zero(3)), constraints(Vector::Zero(3)) {}
    
    Vector& getObjectives() { return objectives; }
    const Vector& getObjectives() const { return objectives; }
    
    Vector& getConstraints() { return constraints; }
    const Vector& getConstraints() const { return constraints; }
    
private:
    Vector objectives;
    Vector constraints;
};

class MockResponseInterface {
public:
    void solve(const MockDesignVariable& design, MockResponse& response) {
        // Simulate solving with simple function
        response.getObjectives() = design.getVector().array().square();
        response.getConstraints() = -design.getVector();
    }
};

class MockDampingInterface {
public:
    void solve(const MockDesignVariable& ref, 
              const Eigen::VectorXd& d,
              const Eigen::VectorXd& rho,
              MockDesignVariable& var,
              MockResponse& response) {
        // Simple approximation for testing
        var = ref;
        response.getObjectives() = ref.getVector().array().square() * 0.9;
        response.getConstraints() = -ref.getVector() * 0.9;
    }
};

class GlobalOptimiserTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize with initial design
        Eigen::Vector3d initialDesign(1.0, 1.0, 1.0);
        designVar = MockDesignVariable(initialDesign);
        
        optimizer = std::make_unique<globopt::GlobalOptimiser<
            MockDesignVariable, MockResponse, MockResponse,
            MockResponseInterface, MockDampingInterface, MockDesignVariable,
            Options>>(designVar, &responseInterface, &dampingInterface);
    }
    
    struct Options {
        Eigen::Vector3d objectiveReference = Eigen::Vector3d::Ones();
        double objectiveTolerence = 1e-6;
        double ConstraintTolerence = 1e-6;
    };
    
    MockDesignVariable designVar;
    MockResponseInterface responseInterface;
    MockDampingInterface dampingInterface;
    std::unique_ptr<globopt::GlobalOptimiser<
        MockDesignVariable, MockResponse, MockResponse,
        MockResponseInterface, MockDampingInterface, MockDesignVariable,
        Options>> optimizer;
        
    const double tolerance = 1e-10;
};

TEST_F(GlobalOptimiserTest, Constructor) {
    EXPECT_TRUE(optimizer->getResponseInterface() == &responseInterface);
    EXPECT_TRUE(optimizer->getDampingInterface() == &dampingInterface);
}

TEST_F(GlobalOptimiserTest, DesignVariableAccess) {
    MockDesignVariable newDesign;
    newDesign.getVector() << 2.0, 2.0, 2.0;
    
    optimizer->setDesignVariable(newDesign);
    auto retrievedDesign = optimizer->getDesignVariable();
    
    EXPECT_TRUE(retrievedDesign.getVector().isApprox(newDesign.getVector(), tolerance));
}

TEST_F(GlobalOptimiserTest, ResponseAccess) {
    // First get initial response
    auto response = optimizer->getResponseProblem();
    
    // Should be zero initially
    EXPECT_TRUE(response.getObjectives().isZero(tolerance));
    EXPECT_TRUE(response.getConstraints().isZero(tolerance));
}

TEST_F(GlobalOptimiserTest, OptimizationFlow) {
    // Setup initial design
    MockDesignVariable initialDesign;
    initialDesign.getVector() << 2.0, 2.0, 2.0;
    optimizer->setDesignVariable(initialDesign);
    
    // Get initial state
    auto initialResponse = optimizer->getResponseProblem();
    
    // Call solve on response interface to update response
    optimizer->getResponseInterface()->solve(initialDesign, initialResponse);
    
    // Run optimization iteration
    MockDesignVariable newDesign;
    newDesign.getVector() << 1.0, 1.0, 1.0;
    optimizer->setDesignVariable(newDesign);
    
    // Call solve again to update response with new design
    MockResponse newResponse;
    optimizer->getResponseInterface()->solve(newDesign, newResponse);
    
    // Check that response changed
    EXPECT_FALSE(newResponse.getObjectives().isApprox(initialResponse.getObjectives()));
    
    // Get final design and check it's different from initial
    auto finalDesign = optimizer->getDesignVariable();
    EXPECT_FALSE(finalDesign.getVector().isApprox(initialDesign.getVector(), tolerance));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
