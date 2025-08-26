/*
 * globalOptimiser.h
 *
 *  Created on: 24 oct. 2016
 *      Author: Ramzi
 */

#ifndef GLOBALOPTIMISER_GLOBALOPTIMISER_HPP_
#define GLOBALOPTIMISER_GLOBALOPTIMISER_HPP_

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <vector>
#include "../GlobalOptimiser/approxFunction.hpp"
#include "dampingUtils.hpp"

namespace globopt {

template<typename designVariable, typename responseProblem, typename dampingProblem,
         typename responseInterface, typename dampingInterface, typename approximationProblem,
		 typename option, typename scalarType = double>
class GlobalOptimiser {
private:
    designVariable designVariables;
    responseProblem responses;
    responseInterface* pResponseInterface;
    responseProblem dampingResponses;
    dampingInterface* pDampingInterface;
    option options;
    
    bool CheckConvergence(const responseProblem& responseCandidate) {
        // Check if max fi < max fref_i for i in O & fi < 0 for i in O_complement
        scalarType relativeObjectiveChange = ((responseCandidate.getObjectives() - responses.getObjectives())
            .cwiseProduct(options.objectiveReference)).maxCoeff();
        scalarType maxConstraint = responseCandidate.getConstraints().maxCoeff();
        
        return (relativeObjectiveChange < options.objectiveTolerence) && (maxConstraint < 0);
    }
    
    bool CheckRelativeConvergence(const responseProblem& responseCandidate) {
        // Check if fi/fref_i - 1 < 0 for i in O
        auto relativeChange = responseCandidate.getObjectives().array() / responses.getObjectives().array() - 1.0;
        return (relativeChange < 0).all();
    }
    
public:
    GlobalOptimiser(designVariable designVariables, 
                   responseInterface* pResponseInterface,  
                   dampingInterface* pDampingInterface) 
        : designVariables(designVariables)
        , pResponseInterface(pResponseInterface)
        , pDampingInterface(pDampingInterface) {}
        
    virtual ~GlobalOptimiser() = default;
    
    void setDesignVariable(const designVariable& variables) {
        this->designVariables = variables;
    }
    
    designVariable getDesignVariable() const {
        return this->designVariables;
    }
    
    responseProblem getResponseProblem() const {
        return this->responses;
    }
    
    responseProblem getResponseDamping() const {
        return this->dampingResponses;
    }
    
    void setResponseInterface(responseInterface* interface) {
        this->pResponseInterface = interface;
    }
    
    responseInterface* getResponseInterface() const {
        return this->pResponseInterface;
    }
    
    void setDampingInterface(dampingInterface* interface) {
        this->pDampingInterface = interface;
    }
    
    dampingInterface* getDampingInterface() const {
        return this->pDampingInterface;
    }
    
    void optimize() {
        // Initial FEM solution
        pResponseInterface->solve(designVariables, responses);
        
        // Store reference point
        designVariable designRef = designVariables;
        responseProblem responseRef = responses;
        
        bool flagOptimiser = false;
        while (!flagOptimiser) {
            // Sub-iteration with approximation
            approximationProblem approx(designRef, options.dampingVector, options.dampingFactors);
            pDampingInterface->solve(approx, designVariables, dampingResponses);
            
            // Update damping factors
            lampar::GlobalOptimiser::UpdateDamping(
                responses.getObjectives(),
                dampingResponses.getObjectives(),
                options.dampingVector[0],
                options.dampingFactors
            );
            
            // FEM solution at new point
            pResponseInterface->solve(designVariables, responses);
            
            // Update damping vector
            lampar::GlobalOptimiser::DampingVector(
                designRef.getVector(),
                designVariables.getVector(),
                options.dampingVector
            );
            
            // Check convergence conditions
            if (CheckConvergence(responses)) {
                designRef = designVariables;
                responseRef = responses;
                
                if (CheckRelativeConvergence(responses)) {
                    flagOptimiser = true;
                }
            }
        }
    }
};

} /* namespace globopt */

#endif /* GLOBALOPTIMISER_GLOBALOPTIMISER_HPP_ */
