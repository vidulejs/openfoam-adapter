#include "Force.H"

using namespace Foam;

preciceAdapter::FSI::Force::Force(
    const Foam::fvMesh& mesh,
    const std::string solverType,
    const std::string nameForce)
: ForceBase(mesh, solverType)
{
    // Check if a force field with the requested name exists.
    // If yes (e.g., solids4Foam), bind Force_ to that field.
    // If not (e.g., pimpleFoam without the Forces function object), create it.
    if (mesh_.foundObject<volVectorField>(nameForce))
    {
        Force_ =
            &const_cast<volVectorField&>(
                mesh_.lookupObject<volVectorField>(nameForce));
    }
    else
    {
        ForceOwning_.reset(new volVectorField(
            IOobject(
                nameForce,
                mesh_.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedVector(
                "fdim",
                dimensionSet(1, 1, -2, 0, 0, 0, 0),
                Foam::vector::zero)));

        Force_ = ForceOwning_.get();
    }
}

// Helper implementation
void preciceAdapter::FSI::Force::resetAccumulator(size_t nData)
{
    accForce_.assign(nData, 0.0);
    accTime_ = 0.0;
}

std::size_t preciceAdapter::FSI::Force::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    std::size_t nData = this->writeToBuffer(buffer, *Force_, dim);

    double currentTime = mesh_.time().value();
    double dt = mesh_.time().deltaTValue();

    // Init buffer if first run
    if (accForce_.size() != nData)
    {
        resetAccumulator(nData);
        lastWriteTime_ = currentTime;

        // Init history vectors
        prevStepForce_.assign(nData, 0.0);
        startWindowForce_.assign(nData, 0.0);
    }
    // 1. Implicit: Time moved backwards
    if (currentTime < lastWriteTime_ - 1e-9)
    {
        // Reset sums
        resetAccumulator(nData);

        prevStepForce_ = startWindowForce_;
    }
    // 2. New Window
    else if (accTime_ >= WINDOW_SIZE - 1e-9)
    {
        startWindowForce_ = prevStepForce_;

        resetAccumulator(nData);
    }

    // Trapezoidal integration
    // Impulse = 1/2 * (F_prev + F_curr) * dt
    for (std::size_t i = 0; i < nData; ++i)
    {
        accForce_[i] += 0.5 * (prevStepForce_[i] + buffer[i]) * dt;

        // Update history for next step
        prevStepForce_[i] = buffer[i];
    }
    accTime_ += dt;

    // OVERWRITE OUTPUT BUFFER WITH AVERAGE
    // F_avg = Sum(Impulse) / Sum(dt)
    if (accTime_ > 1e-12)
    {
        for (std::size_t i = 0; i < nData; ++i)
        {
            buffer[i] = accForce_[i] / accTime_;
        }
    }

    lastWriteTime_ = currentTime;

    return nData;
}

void preciceAdapter::FSI::Force::read(double* buffer, const unsigned int dim)
{
    // Copy the force field from the buffer to OpenFOAM

    // Here we assume that a force volVectorField exists, which is used by
    // the OpenFOAM solver

    int bufferIndex = 0;
    // Set boundary forces
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if (this->locationType_ == LocationType::faceCenters)
        {
            // Make a force field
            vectorField& force = Force_->boundaryFieldRef()[patchID];

            // Copy the forces from the buffer into the force field
            forAll(force, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    force[i][d] = buffer[bufferIndex++];
            }
        }
        else if (this->locationType_ == LocationType::faceNodes)
        {
            // Here we could easily interpolate the face values to point values
            // and assign them to some field, but I guess there is no need
            // unless it will be used
            notImplemented("Read forces not implemented for faceNodes!");
        }
    }

    Foam::scalar maxForceMag = 0;
    Foam::scalar sumForceMag = 0;

    if (this->locationType_ == LocationType::faceCenters)
    {
        for (const label patchID : patchIDs_)
        {
            const fvPatchVectorField& forcePatch = Force_->boundaryField()[patchID];
            forAll(forcePatch, i)
            {
                Foam::scalar magF = mag(forcePatch[i]);
                sumForceMag += magF;
                if (magF > maxForceMag)
                {
                    maxForceMag = magF;
                }
            }
        }
    }

    // write to file
    {
        OFstream os("preciceForceRead.dat", IOstream::ASCII, IOstream::UNCOMPRESSED, IOstream::APPEND);
        os << mesh_.time().timeName() << " "
           << sumForceMag << " "
           << maxForceMag << endl;
    }
}

bool preciceAdapter::FSI::Force::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FSI::Force::getDataName() const
{
    return "Force";
}

Foam::tmp<Foam::vectorField> preciceAdapter::FSI::Force::getFaceVectors(const unsigned int patchID) const
{
    // Normal vectors multiplied by face area
    return mesh_.boundary()[patchID].Sf();
}
