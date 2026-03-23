#include "Laminate/laminateSection.hpp"
#include "Section/section.hpp"

#include <gtest/gtest.h>

#include <memory>

namespace {

template<typename Section>
std::unique_ptr<optsection::section<>> MakeSection(const double lowerThickness,
                                                   const double upperThickness) {
    auto section = std::make_unique<Section>();
    typename Section::Vector_t design = Section::Vector_t::Ones();
    design(Section::Size - 1) = 1.0;
    section->setBoundThickness(lowerThickness, upperThickness);
    section->Initialise(1.0, design);
    return section;
}

}  // namespace

TEST(SectionDispatchTest, DynamicSectionInterfaceHandlesLaminateSections) {
    using UnbalancedSection = lampar::laminateSection<lampar::SingleMaterial, false>;
    using BalancedSection = lampar::laminateSection<lampar::SingleMaterial, true>;

    std::unique_ptr<optsection::section<>> sections[2];
    sections[0] = MakeSection<UnbalancedSection>(0.0, 10.0);
    sections[1] = MakeSection<BalancedSection>(-5.0, 10.0);

    for (auto& section : sections) {
        const int size = section->getSize();
        Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(size, size);
        Eigen::VectorXd design = Eigen::VectorXd::Constant(size, 1.0);
        Eigen::VectorXd residual = Eigen::VectorXd::Ones(size);

        section->HessianEval(hessian);
        EXPECT_EQ(hessian.rows(), size);
        EXPECT_EQ(hessian.cols(), size);
        EXPECT_TRUE(hessian.allFinite());

        const double dualityGap = section->DualityGap();
        EXPECT_GT(dualityGap, 0.0);

        section->CalculateResiduals(design, residual, 0.0);
        EXPECT_EQ(residual.size(), size);
        EXPECT_TRUE(residual.allFinite());

        section->UpdateIncrements(residual);
        double primalStep = 1.0;
        double dualStep = 1.0;
        section->StepSize(primalStep, dualStep);
        EXPECT_GT(primalStep, 0.0);
        EXPECT_GT(dualStep, 0.0);
        EXPECT_LE(primalStep, 1.0);
        EXPECT_LE(dualStep, 1.0);
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
