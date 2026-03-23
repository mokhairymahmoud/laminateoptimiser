#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace lamopt {

using NamedScalarResponseMap = std::unordered_map<std::string, double>;

enum class DerivedScalarResponseTransform {
    Identity,
    Affine,
    AbsoluteAffine,
    InverseAffine
};

struct DerivedScalarResponseRule {
    std::string sourceId;
    DerivedScalarResponseTransform transform = DerivedScalarResponseTransform::Identity;
    double scale = 1.0;
    double offset = 0.0;
    std::string label;
};

[[nodiscard]] inline std::string DescribeDerivedScalarResponseRule(const DerivedScalarResponseRule& rule) {
    if (!rule.label.empty()) {
        return rule.label;
    }
    return rule.sourceId;
}

[[nodiscard]] inline double EvaluateDerivedScalarResponseRule(const DerivedScalarResponseRule& rule,
                                                              const NamedScalarResponseMap& rawValues) {
    const auto iterator = rawValues.find(rule.sourceId);
    if (iterator == rawValues.end()) {
        throw std::runtime_error("Response schema source id was not extracted: " + rule.sourceId);
    }

    const double value = iterator->second;
    switch (rule.transform) {
        case DerivedScalarResponseTransform::Identity:
            return value;
        case DerivedScalarResponseTransform::Affine:
            return rule.scale * value + rule.offset;
        case DerivedScalarResponseTransform::AbsoluteAffine:
            return rule.scale * std::abs(value) + rule.offset;
        case DerivedScalarResponseTransform::InverseAffine:
            if (std::abs(value) <= 1.0e-16) {
                throw std::runtime_error("Response schema inverse transform encountered a zero source value for "
                                         + DescribeDerivedScalarResponseRule(rule) + ".");
            }
            return rule.scale / value + rule.offset;
    }

    throw std::runtime_error("Unsupported response schema transform for "
                             + DescribeDerivedScalarResponseRule(rule) + ".");
}

[[nodiscard]] inline double EvaluateDerivedScalarResponseDerivative(const DerivedScalarResponseRule& rule,
                                                                    const NamedScalarResponseMap& rawValues) {
    const auto iterator = rawValues.find(rule.sourceId);
    if (iterator == rawValues.end()) {
        throw std::runtime_error("Response schema source id was not extracted: " + rule.sourceId);
    }

    const double value = iterator->second;
    switch (rule.transform) {
        case DerivedScalarResponseTransform::Identity:
            return 1.0;
        case DerivedScalarResponseTransform::Affine:
            return rule.scale;
        case DerivedScalarResponseTransform::AbsoluteAffine:
            if (value > 0.0) {
                return rule.scale;
            }
            if (value < 0.0) {
                return -rule.scale;
            }
            return 0.0;
        case DerivedScalarResponseTransform::InverseAffine:
            if (std::abs(value) <= 1.0e-16) {
                throw std::runtime_error("Response schema inverse derivative encountered a zero source value for "
                                         + DescribeDerivedScalarResponseRule(rule) + ".");
            }
            return -rule.scale / (value * value);
    }

    throw std::runtime_error("Unsupported response schema transform derivative for "
                             + DescribeDerivedScalarResponseRule(rule) + ".");
}

}  // namespace lamopt
