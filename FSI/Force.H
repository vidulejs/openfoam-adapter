#ifndef FSI_FORCE_H
#define FSI_FORCE_H

#include "ForceBase.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads force
class Force : public ForceBase
{
private:
    //- Force field (non-owning pointer, it may already be constructed by the solver)
    Foam::volVectorField* Force_;

    //- Force field (owning pointer, we bind to Force_)
    std::unique_ptr<Foam::volVectorField> ForceOwning_;

    std::vector<double> accForce_; // Sum of Force * dt
    double accTime_ = 0.0;         // Sum of dt
    double lastWriteTime_ = -1.0;

    void resetAccumulator(size_t nData);
    // HARDCODED Window Size
    // this needs to be read programmatically from another class
    const double WINDOW_SIZE = 0.01;
    std::vector<double> prevStepForce_;    // prev substep force
    std::vector<double> startWindowForce_; // Force at Window Start

public:
    //- Constructor
    Force(
        const Foam::fvMesh& mesh,
        const std::string solverType,
        const std::string nameForce);

    //- Write the forces values into the buffer
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the forces values from the buffer
    void read(double* buffer, const unsigned int dim) final;

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- Returns the normal vectors multiplied by the face area
    Foam::tmp<Foam::vectorField> getFaceVectors(const unsigned int patchID) const final;
};

}
}

#endif
