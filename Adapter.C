#include "Adapter.H"
#include "Interface.H"
#include "Utilities.H"

#include "IOstreams.H"

using namespace Foam;

preciceAdapter::Adapter::Adapter(const Time& runTime, const fvMesh& mesh)
: runTime_(runTime),
  mesh_(mesh)
{
    adapterInfo("Loaded the OpenFOAM-preCICE adapter - v1.3.1.", "info");

    return;
}

bool preciceAdapter::Adapter::configFileRead()
{

    // We need a try-catch here, as if reading preciceDict fails,
    // the respective exception will be reduced to a warning.
    // See also comment in preciceAdapter::Adapter::configure().
    try
    {
        SETUP_TIMER();
        adapterInfo("Reading preciceDict...", "info");

        // TODO: static is just a quick workaround to be able
        // to find the dictionary also out of scope (e.g. in KappaEffective).
        // We need a better solution.
        static IOdictionary preciceDict(
            IOobject(
                "preciceDict",
                runTime_.system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE));

        // Read and display the preCICE configuration file name
        preciceConfigFilename_ = preciceDict.get<fileName>("preciceConfig");
        DEBUG(adapterInfo("  precice-config-file : " + preciceConfigFilename_));

        // Read and display the participant name
        participantName_ = preciceDict.get<word>("participant");
        DEBUG(adapterInfo("  participant name    : " + participantName_));

        // Read and display the list of modules
        DEBUG(adapterInfo("  modules requested   : "));
        auto modules_ = preciceDict.get<wordList>("modules");
        for (const auto& module : modules_)
        {
            DEBUG(adapterInfo("  - " + module + "\n"));

            // Set the modules switches
            if (module == "CHT")
            {
                CHTenabled_ = true;
            }

            if (module == "FSI")
            {
                FSIenabled_ = true;
            }

            if (module == "FF")
            {
                FFenabled_ = true;
            }
        }

        // Every interface is a subdictionary of "interfaces",
        // each with an arbitrary name. Read all of them and create
        // a list (here: pointer) of dictionaries.
        const auto* interfaceDictPtr = preciceDict.findDict("interfaces");
        DEBUG(adapterInfo("  interfaces : "));

        // Check if we found any interfaces
        // and get the details of each interface
        if (!interfaceDictPtr)
        {
            adapterInfo("  Empty list of interfaces", "warning");
            return false;
        }
        else
        {
            for (const entry& interfaceDictEntry : *interfaceDictPtr)
            {
                if (interfaceDictEntry.isDict())
                {
                    const dictionary& interfaceDict = interfaceDictEntry.dict();
                    struct InterfaceConfig interfaceConfig;

                    interfaceConfig.meshName = interfaceDict.get<word>("mesh");
                    DEBUG(adapterInfo("  - mesh         : " + interfaceConfig.meshName));

                    // By default, assume "faceCenters" as locationsType
                    interfaceConfig.locationsType = interfaceDict.lookupOrDefault<word>("locations", "faceCenters");
                    DEBUG(adapterInfo("    locations    : " + interfaceConfig.locationsType));

                    // By default, assume that no mesh connectivity is required (i.e. no nearest-projection mapping)
                    interfaceConfig.meshConnectivity = interfaceDict.lookupOrDefault<bool>("connectivity", false);
                    // Mesh connectivity only makes sense in case of faceNodes, check and raise a warning otherwise
                    if (interfaceConfig.meshConnectivity && (interfaceConfig.locationsType == "faceCenters" || interfaceConfig.locationsType == "volumeCenters" || interfaceConfig.locationsType == "volumeCentres"))
                    {
                        DEBUG(adapterInfo("Mesh connectivity is not supported for faceCenters or volumeCenters. \n"
                                          "Please configure the desired interface with the locationsType faceNodes. \n"
                                          "Have a look in the adapter documentation for detailed information.",
                                          "warning"));
                        return false;
                    }
                    DEBUG(adapterInfo("    connectivity : " + std::to_string(interfaceConfig.meshConnectivity)));

                    DEBUG(adapterInfo("    patches      : "));
                    auto patches = interfaceDict.get<wordList>("patches");
                    for (auto patch : patches)
                    {
                        interfaceConfig.patchNames.push_back(patch);
                        DEBUG(adapterInfo("      - " + patch));
                    }

                    DEBUG(adapterInfo("    cellSets      : "));
                    auto cellSets = interfaceDict.lookupOrDefault<wordList>("cellSets", wordList());

                    for (auto cellSet : cellSets)
                    {
                        interfaceConfig.cellSetNames.push_back(cellSet);
                        DEBUG(adapterInfo("      - " + cellSet));
                    }

                    if (!interfaceConfig.cellSetNames.empty() && !(interfaceConfig.locationsType == "volumeCenters" || interfaceConfig.locationsType == "volumeCentres"))
                    {
                        adapterInfo("Cell sets are not supported for locationType != volumeCenters. \n"
                                    "Please configure the desired interface with the locationsType volumeCenters. \n"
                                    "Have a look in the adapter documentation for detailed information.",
                                    "warning");
                        return false;
                    }

                    DEBUG(adapterInfo("    writeData    : "));
                    auto writeData = interfaceDict.get<wordList>("writeData");
                    for (auto writeDatum : writeData)
                    {
                        interfaceConfig.writeData.push_back(writeDatum);
                        DEBUG(adapterInfo("      - " + writeDatum));
                    }

                    DEBUG(adapterInfo("    readData     : "));
                    auto readData = interfaceDict.get<wordList>("readData");
                    for (auto readDatum : readData)
                    {
                        interfaceConfig.readData.push_back(readDatum);
                        DEBUG(adapterInfo("      - " + readDatum));
                    }
                    interfacesConfig_.push_back(interfaceConfig);
                }
            }
        }

        // NOTE: set the switch for your new module here

        // If the CHT module is enabled, create it, read the
        // CHT-specific options and configure it.
        if (CHTenabled_)
        {
            CHT_ = new CHT::ConjugateHeatTransfer(mesh_);
            if (!CHT_->configure(preciceDict))
            {
                return false;
            }
        }

        // If the FSI module is enabled, create it, read the
        // FSI-specific options and configure it.
        if (FSIenabled_)
        {
            // Check for unsupported FSI with meshConnectivity
            for (uint i = 0; i < interfacesConfig_.size(); i++)
            {
                if (interfacesConfig_.at(i).meshConnectivity == true)
                {
                    adapterInfo(
                        "You have requested mesh connectivity (most probably for nearest-projection mapping) "
                        "and you have enabled the FSI module. "
                        "Mapping with connectivity information is not implemented for FSI, only for CHT-related fields. "
                        "warning");
                    return false;
                }
            }

            FSI_ = new FSI::FluidStructureInteraction(mesh_, runTime_);
            if (!FSI_->configure(preciceDict))
            {
                return false;
            }
        }

        if (FFenabled_)
        {
            FF_ = new FF::FluidFluid(mesh_);
            if (!FF_->configure(preciceDict))
            {
                return false;
            }
        }

        // NOTE: Create your module and read any options specific to it here

        if (!CHTenabled_ && !FSIenabled_ && !FFenabled_) // NOTE: Add your new switch here
        {
            adapterInfo("No module is enabled.", "error-deferred");
            return false;
        }

        // TODO: Loading modules should be implemented in more general way,
        // in order to avoid code duplication. See issue #16 on GitHub.

        ACCUMULATE_TIMER(timeInConfigRead_);
    }
    catch (const Foam::error& e)
    {
        adapterInfo(e.message(), "error-deferred");
        return false;
    }

    return true;
}

void preciceAdapter::Adapter::configure()
{
    // Read the adapter's configuration file
    if (!configFileRead())
    {
        // This method is called from the functionObject's read() method,
        // which is called by the Foam::functionObjectList::read() method.
        // All the exceptions triggered in this method are caught as
        // warnings and the simulation continues simply without the
        // functionObject. However, we want the simulation to exit with an
        // error in case something is wrong. We store the information that
        // there was an error and it will be handled by the first call to
        // the functionObject's execute(), which can throw errors normally.
        errorsInConfigure = true;

        return;
    }

    try
    {
        // Check the timestep type (fixed vs adjustable)
        DEBUG(adapterInfo("Checking the timestep type (fixed vs adjustable)..."));
        adjustableTimestep_ = runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);

        if (adjustableTimestep_)
        {
            DEBUG(adapterInfo("  Timestep type: adjustable."));
        }
        else
        {
            DEBUG(adapterInfo("  Timestep type: fixed."));
        }

        // Construct preCICE
        SETUP_TIMER();
        DEBUG(adapterInfo("Creating the preCICE solver interface..."));
        DEBUG(adapterInfo("  Number of processes: " + std::to_string(Pstream::nProcs())));
        DEBUG(adapterInfo("  MPI rank: " + std::to_string(Pstream::myProcNo())));
        precice_ = new precice::Participant(participantName_, preciceConfigFilename_, Pstream::myProcNo(), Pstream::nProcs());
        DEBUG(adapterInfo("  preCICE solver interface was created."));

        ACCUMULATE_TIMER(timeInPreciceConstruct_);

        // Create interfaces
        REUSE_TIMER();
        DEBUG(adapterInfo("Creating interfaces..."));
        for (uint i = 0; i < interfacesConfig_.size(); i++)
        {
            std::string namePointDisplacement = FSIenabled_ ? FSI_->getPointDisplacementFieldName() : "default";
            std::string nameCellDisplacement = FSIenabled_ ? FSI_->getCellDisplacementFieldName() : "default";
            bool restartFromDeformed = FSIenabled_ ? FSI_->isRestartingFromDeformed() : false;

            Interface* interface = new Interface(*precice_, mesh_, interfacesConfig_.at(i).meshName, interfacesConfig_.at(i).locationsType, interfacesConfig_.at(i).patchNames, interfacesConfig_.at(i).cellSetNames, interfacesConfig_.at(i).meshConnectivity, restartFromDeformed, namePointDisplacement, nameCellDisplacement);
            interfaces_.push_back(interface);
            DEBUG(adapterInfo("Interface created on mesh " + interfacesConfig_.at(i).meshName));

            DEBUG(adapterInfo("Adding coupling data writers..."));
            for (uint j = 0; j < interfacesConfig_.at(i).writeData.size(); j++)
            {
                std::string dataName = interfacesConfig_.at(i).writeData.at(j);

                unsigned int inModules = 0;

                // Add CHT-related coupling data writers
                if (CHTenabled_ && CHT_->addWriters(dataName, interface))
                {
                    inModules++;
                }

                // Add FSI-related coupling data writers
                if (FSIenabled_ && FSI_->addWriters(dataName, interface))
                {
                    inModules++;
                }

                // Add FF-related coupling data writers
                if (FFenabled_ && FF_->addWriters(dataName, interface))
                {
                    inModules++;
                }

                if (inModules == 0)
                {
                    adapterInfo("I don't know how to write \"" + dataName
                                    + "\". Maybe this is a typo or maybe you need to enable some adapter module?",
                                "error-deferred");
                }
                else if (inModules > 1)
                {
                    adapterInfo("It looks like more than one modules can write \"" + dataName
                                    + "\" and I don't know how to choose. Try disabling one of the modules.",
                                "error-deferred");
                }

                // NOTE: Add any coupling data writers for your module here.
            } // end add coupling data writers

            DEBUG(adapterInfo("Adding coupling data readers..."));
            for (uint j = 0; j < interfacesConfig_.at(i).readData.size(); j++)
            {
                std::string dataName = interfacesConfig_.at(i).readData.at(j);

                unsigned int inModules = 0;

                // Add CHT-related coupling data readers
                if (CHTenabled_ && CHT_->addReaders(dataName, interface)) inModules++;

                // Add FSI-related coupling data readers
                if (FSIenabled_ && FSI_->addReaders(dataName, interface)) inModules++;

                // Add FF-related coupling data readers
                if (FFenabled_ && FF_->addReaders(dataName, interface)) inModules++;

                if (inModules == 0)
                {
                    adapterInfo("I don't know how to read \"" + dataName
                                    + "\". Maybe this is a typo or maybe you need to enable some adapter module?",
                                "error-deferred");
                }
                else if (inModules > 1)
                {
                    adapterInfo("It looks like more than one modules can read \"" + dataName
                                    + "\" and I don't know how to choose. Try disabling one of the modules.",
                                "error-deferred");
                }

                // NOTE: Add any coupling data readers for your module here.
            } // end add coupling data readers

            // Create the interface's data buffer
            interface->createBuffer();
        }
        ACCUMULATE_TIMER(timeInMeshSetup_);

        // Initialize preCICE and exchange the first coupling data
        initialize();

        // If checkpointing is required, specify the checkpointed fields
        // and write the first checkpoint
        if (requiresWritingCheckpoint())
        {
            checkpointing_ = true;

            // Setup the checkpointing (find and add fields to checkpoint)
            setupCheckpointing();

            // Write checkpoint (for the first iteration)
            writeCheckpoint();
        }

        // Adjust the timestep for the first iteration, if it is fixed
        if (!adjustableTimestep_)
        {
            adjustSolverTimeStepAndReadData();
        }

        // If the solver tries to end before the coupling is complete,
        // e.g. because the solver's endTime was smaller or (in implicit
        // coupling) equal with the max-time specified in preCICE,
        // problems may occur near the end of the simulation,
        // as the function object may be called only once near the end.
        // See the implementation of Foam::Time::run() for more details.
        // To prevent this, we set the solver's endTime to "infinity"
        // and let only preCICE control the end of the simulation.
        // This has the side-effect of not triggering the end() method
        // in any function object normally. Therefore, we trigger it
        // when preCICE dictates to stop the coupling.
        adapterInfo(
            "Setting the solver's endTime to infinity to prevent early exits. "
            "Only preCICE will control the simulation's endTime. "
            "Any functionObject's end() method will be triggered by the adapter. "
            "You may disable this behavior in the adapter's configuration.",
            "info");
        const_cast<Time&>(runTime_).setEndTime(GREAT);
    }
    catch (const Foam::error& e)
    {
        adapterInfo(e.message(), "error-deferred");
        errorsInConfigure = true;
    }

    return;
}

void preciceAdapter::Adapter::execute()
{
    if (errorsInConfigure)
    {
        // Handle any errors during configure().
        // See the comments in configure() for details.
        adapterInfo(
            "There was a problem while configuring the adapter. "
            "See the log for details.",
            "error");
    }

    // The solver has already solved the equations for this timestep.
    // Now call the adapter's methods to perform the coupling.

    // TODO add a function which checks if all fields are checkpointed.
    // if (ncheckpointed is nregisterdobjects. )

    // Write the coupling data in the buffer
    writeCouplingData();

    // Advance preCICE
    advance();

    // Read checkpoint if required
    if (requiresReadingCheckpoint())
    {
        pruneCheckpointedFields();
        readCheckpoint();
    }

    // Write checkpoint if required
    if (requiresWritingCheckpoint())
    {
        writeCheckpoint();
    }

    // As soon as OpenFOAM writes the results, it will not try to write again
    // if the time takes the same value again. Therefore, during an implicit
    // coupling, we write again when the coupling timestep is complete.
    // Check the behavior e.g. by using watch on a result file:
    //     watch -n 0.1 -d ls --full-time Fluid/0.01/T.gz
    SETUP_TIMER();
    if (checkpointing_ && isCouplingTimeWindowComplete())
    {
        // Check if the time directory already exists
        // (i.e. the solver wrote results that need to be updated)
        if (runTime_.timePath().type() == fileName::DIRECTORY)
        {
            adapterInfo(
                "The coupling timestep completed. "
                "Writing the updated results.",
                "info");
            const_cast<Time&>(runTime_).writeNow();
        }
    }
    ACCUMULATE_TIMER(timeInWriteResults_);

    // Adjust the timestep, if it is fixed
    if (!adjustableTimestep_)
    {
        adjustSolverTimeStepAndReadData();
    }

    // If the coupling is not going to continue, tear down everything
    // and stop the simulation.
    if (!isCouplingOngoing())
    {
        adapterInfo("The coupling completed.", "info");

        // Finalize the preCICE solver interface and delete data
        finalize();

        // Tell OpenFOAM to stop the simulation.
        // Set the solver's endTime to now. The next evaluation of
        // runTime.run() will be false and the solver will exit.
        const_cast<Time&>(runTime_).setEndTime(runTime_.value());
        adapterInfo(
            "The simulation was ended by preCICE. "
            "Calling the end() methods of any functionObject explicitly.",
            "info");
        adapterInfo("Great that you are using the OpenFOAM-preCICE adapter! "
                    "Next to the preCICE library and any other components, please also cite this adapter. "
                    "Find how on https://precice.org/adapter-openfoam-overview.html.",
                    "info");
        const_cast<Time&>(runTime_).functionObjects().end();
    }

    return;
}


void preciceAdapter::Adapter::adjustTimeStep()
{
    adjustSolverTimeStepAndReadData();

    return;
}

void preciceAdapter::Adapter::readCouplingData(double relativeReadTime)
{
    SETUP_TIMER();
    DEBUG(adapterInfo("Reading coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->readCouplingData(relativeReadTime);
    }

    ACCUMULATE_TIMER(timeInRead_);

    return;
}

void preciceAdapter::Adapter::writeCouplingData()
{
    SETUP_TIMER();
    DEBUG(adapterInfo("Writing coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->writeCouplingData();
    }

    ACCUMULATE_TIMER(timeInWrite_);

    return;
}

void preciceAdapter::Adapter::initialize()
{
    DEBUG(adapterInfo("Initializing the preCICE solver interface..."));
    SETUP_TIMER();

    if (precice_->requiresInitialData())
        writeCouplingData();

    DEBUG(adapterInfo("Initializing preCICE data..."));
    precice_->initialize();
    preciceInitialized_ = true;
    ACCUMULATE_TIMER(timeInInitialize_);

    adapterInfo("preCICE was configured and initialized", "info");

    return;
}

void preciceAdapter::Adapter::finalize()
{
    if (NULL != precice_ && preciceInitialized_ && !isCouplingOngoing())
    {
        DEBUG(adapterInfo("Finalizing the preCICE solver interface..."));

        // Finalize the preCICE solver interface
        SETUP_TIMER();
        precice_->finalize();
        ACCUMULATE_TIMER(timeInFinalize_);

        preciceInitialized_ = false;

        // Delete the solver interface and all the related data
        teardown();
    }
    else
    {
        adapterInfo("Could not finalize preCICE.", "error");
    }

    return;
}

void preciceAdapter::Adapter::advance()
{
    DEBUG(adapterInfo("Advancing preCICE..."));

    SETUP_TIMER();
    precice_->advance(timestepSolver_);
    ACCUMULATE_TIMER(timeInAdvance_);

    return;
}

void preciceAdapter::Adapter::adjustSolverTimeStepAndReadData()
{
    DEBUG(adapterInfo("Adjusting the solver's timestep..."));

    // The timestep size that the solver has determined that it wants to use
    double timestepSolverDetermined;

    /* In this method, the adapter overwrites the timestep used by OpenFOAM.
       If the timestep is not adjustable, OpenFOAM will not try to re-estimate
       the timestep or read it again from the controlDict. Therefore, store
       the value that the timestep has is the beginning and try again to use this
       in every iteration.
       // TODO Treat also the case where the user modifies the timestep
       // in the controlDict during the simulation.
    */

    // Is the timestep adjustable or fixed?
    if (!adjustableTimestep_)
    {
        // Have we already stored the timestep?
        if (!useStoredTimestep_)
        {
            // Show a warning if runTimeModifiable is set
            if (runTime_.runTimeModifiable())
            {
                adapterInfo(
                    "You have enabled 'runTimeModifiable' in the "
                    "controlDict. The preciceAdapter does not yet "
                    "fully support this functionality when "
                    "'adjustableTimestep' is not enabled. "
                    "If you modify the 'deltaT' in the controlDict "
                    "during the simulation, it will not be updated.",
                    "warning");
            }

            // Store the value
            timestepStored_ = runTime_.deltaT().value();

            // Ok, we stored it once, we will use this from now on
            useStoredTimestep_ = true;
        }

        // Use the stored timestep as the determined solver's timestep
        timestepSolverDetermined = timestepStored_;
    }
    else
    {
        // The timestep is adjustable, so OpenFOAM will modify it
        // and therefore we can use the updated value
        timestepSolverDetermined = runTime_.deltaT().value();
    }

    /* If the solver tries to use a timestep smaller than the one determined
       by preCICE, that means that the solver is trying to subcycle.
       This may not be allowed by the user.
       If the solver tries to use a bigger timestep, then it needs to use
       the same timestep as the one determined by preCICE.
    */
    double tolerance = 1e-14;
    if (precice_->getMaxTimeStepSize() - timestepSolverDetermined > tolerance)
    {
        // Add a bool 'subCycling = true' which is checked in the storeMeshPoints() function.
        adapterInfo(
            "The solver's timestep is smaller than the "
            "coupling timestep. Subcycling...",
            "info");
        timestepSolver_ = timestepSolverDetermined;
        // TODO subcycling is enabled. For FSI the oldVolumes must be written, which is normally not done.
        if (FSIenabled_)
        {
            adapterInfo(
                "The adapter does not fully support subcycling for FSI and instabilities may occur.",
                "warning");
        }
    }
    else if (timestepSolverDetermined - precice_->getMaxTimeStepSize() > tolerance)
    {
        // In the last time-step, we adjust to dt = 0, but we don't need to trigger the warning here
        if (precice_->isCouplingOngoing())
        {
            adapterInfo(
                "The solver's timestep cannot be larger than the coupling timestep."
                " Adjusting from "
                    + std::to_string(timestepSolverDetermined) + " to " + std::to_string(precice_->getMaxTimeStepSize()),
                "warning");
        }
        timestepSolver_ = precice_->getMaxTimeStepSize();
    }
    else
    {
        DEBUG(adapterInfo("The solver's timestep is the same as the "
                          "coupling timestep."));
        timestepSolver_ = precice_->getMaxTimeStepSize();
    }

    // Update the solver's timestep (but don't trigger the adjustDeltaT(),
    // which also triggers the functionObject's adjustTimeStep())
    // TODO: Keep this in mind if any relevant problem appears.
    const_cast<Time&>(runTime_).setDeltaT(timestepSolver_, false);

    DEBUG(adapterInfo("Reading coupling data associated to the calculated time-step size..."));

    // Read the received coupling data from the buffer
    // Fits to an implicit Euler
    readCouplingData(runTime_.deltaT().value());

    return;
}

bool preciceAdapter::Adapter::isCouplingOngoing()
{
    bool isCouplingOngoing = false;

    // If the coupling ends before the solver ends,
    // the solver would try to access this method again,
    // giving a segmentation fault if precice_
    // was not available.
    if (NULL != precice_)
    {
        isCouplingOngoing = precice_->isCouplingOngoing();
    }

    return isCouplingOngoing;
}

bool preciceAdapter::Adapter::isCouplingTimeWindowComplete()
{
    return precice_->isTimeWindowComplete();
}

bool preciceAdapter::Adapter::requiresReadingCheckpoint()
{
    return precice_->requiresReadingCheckpoint();
}

bool preciceAdapter::Adapter::requiresWritingCheckpoint()
{
    return precice_->requiresWritingCheckpoint();
}


void preciceAdapter::Adapter::storeCheckpointTime()
{
    couplingIterationTimeIndex_ = runTime_.timeIndex();
    couplingIterationTimeValue_ = runTime_.value();
    DEBUG(adapterInfo("Stored time value t = " + std::to_string(runTime_.value())));

    return;
}

void preciceAdapter::Adapter::reloadCheckpointTime()
{
    const_cast<Time&>(runTime_).setTime(couplingIterationTimeValue_, couplingIterationTimeIndex_);
    // TODO also reset the current iteration?!
    DEBUG(adapterInfo("Reloaded time value t = " + std::to_string(runTime_.value())));

    return;
}

void preciceAdapter::Adapter::storeMeshPoints()
{
    DEBUG(adapterInfo("Storing mesh points..."));
    // TODO: In foam-extend, we would need "allPoints()". Check if this gives the same data.
    meshPoints_ = mesh_.points();
    oldMeshPoints_ = mesh_.oldPoints();

    /*
    // TODO  This is only required for subcycling. It should not be called when not subcycling!!
    // Add a bool 'subcycling' which can be evaluated every timestep.
    if ( !oldVolsStored && mesh_.foundObject<volScalarField::Internal>("V00") ) // For Ddt schemes which use one previous timestep
    {
        setupMeshVolCheckpointing();
        oldVolsStored = true;
    }
    // Update any volume fields from the buffer to the checkpointed values (if already exists.)
    */

    DEBUG(adapterInfo("Stored mesh points."));
    if (mesh_.moving())
    {
        if (!meshCheckPointed)
        {
            // Set up the checkpoint for the mesh flux: meshPhi
            setupMeshCheckpointing();
            meshCheckPointed = true;
        }
        writeMeshCheckpoint();
        writeVolCheckpoint(); // Does not write anything unless subcycling.
    }
}

void preciceAdapter::Adapter::reloadMeshPoints()
{
    if (!mesh_.moving())
    {
        DEBUG(adapterInfo("Mesh points not moved as the mesh is not moving"));
        return;
    }

    // In Foam::polyMesh::movePoints.
    // TODO: The function movePoints overwrites the pointer to the old mesh.
    // Therefore, if you revert the mesh, the oldpointer will be set to the points, which are the new values.
    DEBUG(adapterInfo("Moving mesh points to their previous locations..."));

    // TODO
    // Switch oldpoints on for pure physics. (is this required?). Switch off for better mesh deformation capabilities?
    // const_cast<pointField&>(mesh_.points()) = oldMeshPoints_;
    const_cast<fvMesh&>(mesh_).movePoints(meshPoints_);

    DEBUG(adapterInfo("Moved mesh points to their previous locations."));

    // TODO The if statement can be removed in this case, but it is still included for clarity
    if (meshCheckPointed)
    {
        readMeshCheckpoint();
    }

    /*  // TODO This part should only be used when sybcycling. See the description in 'storeMeshPoints()'
        // The if statement can be removed in this case, but it is still included for clarity
    if ( oldVolsStored )
    {
        readVolCheckpoint();
    }
    */
}

void preciceAdapter::Adapter::setupMeshCheckpointing()
{
    // The other mesh <type>Fields:
    //      C
    //      Cf
    //      Sf
    //      magSf
    //      delta
    // are updated by the function fvMesh::movePoints. Only the meshPhi needs checkpointing.
    DEBUG(adapterInfo("Creating a list of the mesh checkpointed fields..."));

    // Add meshPhi to the checkpointed fields
    addMeshCheckpointField(
        const_cast<surfaceScalarField&>(
            mesh_.phi()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.phi().name() + " in the list of checkpointed fields.");
#endif
}

void preciceAdapter::Adapter::setupMeshVolCheckpointing()
{
    DEBUG(adapterInfo("Creating a list of the mesh volume checkpointed fields..."));
    // Add the V0 and the V00 to the list of checkpointed fields.
    // For V0
    addVolCheckpointField(
        const_cast<volScalarField::Internal&>(
            mesh_.V0()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of checkpointed fields.");
#endif
    // For V00
    addVolCheckpointField(
        const_cast<volScalarField::Internal&>(
            mesh_.V00()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of checkpointed fields.");
#endif

    // Also add the buffer fields.
    // TODO For V0
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V0()
        )
    ); */
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of buffer checkpointed fields.");
#endif
    // TODO For V00
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V00()
        )
    );*/
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of buffer checkpointed fields.");
#endif
}


void preciceAdapter::Adapter::setupCheckpointing()
{
    SETUP_TIMER();

    // Add fields in the checkpointing list - sorted for parallel consistency
    DEBUG(adapterInfo("Adding in checkpointed fields..."));

#undef doLocalCode
#define doLocalCode(GeomField)                                           \
    /* Checkpoint registered GeomField objects */                        \
    for (const word& obj : mesh_.sortedNames<GeomField>())               \
    {                                                                    \
        std::string objStr = obj;                                        \
        DEBUG(adapterInfo("Checkpoint " + objStr + " : " #GeomField));      \
        addCheckpointField(objStr, mesh_.thisDb().getObjectPtr<GeomField>(obj)); \
    }

    doLocalCode(volScalarField);
    doLocalCode(volVectorField);
    doLocalCode(volTensorField);
    doLocalCode(volSymmTensorField);

    doLocalCode(surfaceScalarField);
    doLocalCode(surfaceVectorField);
    doLocalCode(surfaceTensorField);

    doLocalCode(pointScalarField);
    doLocalCode(pointVectorField);
    doLocalCode(pointTensorField);

    // NOTE: Add here other object types to checkpoint, if needed.

#undef doLocalCode

    ACCUMULATE_TIMER(timeInCheckpointingSetup_);
}

void preciceAdapter::Adapter::printFieldCountsDB()
{

// Print fields in the OpenFOAM database
// don't print the fields ending with '0' aka oldTime

std::string fields;
int count = 0;
int strlen = 0;

#undef doLocalCode
#define doLocalCode(GeomField)                                           \
    for (const word& obj : mesh_.sortedNames<GeomField>())               \
    {                                                                    \
        strlen = obj.size();                                             \
        if (obj[strlen - 1] != '0')                                      \
        {                                                                \
            fields += obj + " ";                                         \
            count++;                                                     \
        }                                                                \
    }                                                                    \
    DEBUG(adapterInfo(std::to_string(count) + " : " #GeomField + " : " + fields));
    fields = "";
    count = 0;

    doLocalCode(volScalarField);
    doLocalCode(volVectorField);
    doLocalCode(volTensorField);
    doLocalCode(volSymmTensorField);

    doLocalCode(surfaceScalarField);
    doLocalCode(surfaceVectorField);
    doLocalCode(surfaceTensorField);

    doLocalCode(pointScalarField);
    doLocalCode(pointVectorField);
    doLocalCode(pointTensorField);

#undef doLocalCode
}

void preciceAdapter::Adapter::pruneCheckpointedFields()
{
    // Check if checkpointed fields exist in OpenFOAM database
    // If not, remove them from the checkpointed fields list

    word obj;
    std::vector<word> fieldNames;
    std::vector<std::string> toRemove;
#undef doLocalCode
#define doLocalCode(GeomField, GeomFieldMap, GeomFieldCopiesMap)           \
    fieldNames.clear();                                                    \
    toRemove.clear();                                                      \
    for (const word& obj : mesh_.sortedNames<GeomField>())                 \
    {                                                                      \
        fieldNames.push_back(obj);                                         \
    }                                                                      \
    for (const auto& kv : GeomFieldMap){                                   \
        obj = kv.first;                                                    \
        if (std::find(fieldNames.begin(), fieldNames.end(), obj) == fieldNames.end())  \
        {                                                                  \
            toRemove.push_back(static_cast<std::string>(obj));             \
        }                                                                  \
    }                                                                      \
    for (const auto& obj : toRemove) {                                     \
        GeomFieldMap.erase(obj);                                           \
        delete GeomFieldCopiesMap[obj];                                    \
        GeomFieldCopiesMap.erase(obj);                                     \
        DEBUG(adapterInfo("Removed " #GeomField " : " + obj + " from the checkpointed fields list.")); \
    }

doLocalCode(volScalarField, volScalarFields_, volScalarFieldCopies_);
doLocalCode(volVectorField, volVectorFields_, volVectorFieldCopies_);
doLocalCode(volTensorField, volTensorFields_, volTensorFieldCopies_);
doLocalCode(volSymmTensorField, volSymmTensorFields_, volSymmTensorFieldCopies_);

doLocalCode(surfaceScalarField, surfaceScalarFields_, surfaceScalarFieldCopies_);
doLocalCode(surfaceVectorField, surfaceVectorFields_, surfaceVectorFieldCopies_);
doLocalCode(surfaceTensorField, surfaceTensorFields_, surfaceTensorFieldCopies_);

doLocalCode(pointScalarField, pointScalarFields_, pointScalarFieldCopies_);
doLocalCode(pointVectorField, pointVectorFields_, pointVectorFieldCopies_);
doLocalCode(pointTensorField, pointTensorFields_, pointTensorFieldCopies_);

#undef doLocalCode
}

// All mesh checkpointed fields

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceScalarField& field)
{
    {
        meshSurfaceScalarFields_.push_back(&field);
        meshSurfaceScalarFieldCopies_.push_back(new surfaceScalarField(field));
    }
}

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceVectorField& field)
{
    {
        meshSurfaceVectorFields_.push_back(&field);
        meshSurfaceVectorFieldCopies_.push_back(new surfaceVectorField(field));
    }
}

void preciceAdapter::Adapter::addMeshCheckpointField(volVectorField& field)
{
    {
        meshVolVectorFields_.push_back(&field);
        meshVolVectorFieldCopies_.push_back(new volVectorField(field));
    }
}

// TODO Internal field for the V0 (volume old) and V00 (volume old-old) fields
void preciceAdapter::Adapter::addVolCheckpointField(volScalarField::Internal& field)
{
    {
        volScalarInternalFields_.push_back(&field);
        volScalarInternalFieldCopies_.push_back(new volScalarField::Internal(field));
    }
}


void preciceAdapter::Adapter::addCheckpointField(const std::string& name, volScalarField* field)
{
    if (field)
    {
        volScalarFields_[name] = field;
        volScalarFieldCopies_[name] = new volScalarField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, volVectorField* field)
{
    if (field)
    {
        volVectorFields_[name] = field;
        volVectorFieldCopies_[name] = new volVectorField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, surfaceScalarField* field)
{
    if (field)
    {
        surfaceScalarFields_[name] = field;
        surfaceScalarFieldCopies_[name] = new surfaceScalarField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, surfaceVectorField* field)
{
    if (field)
    {
        surfaceVectorFields_[name] = field;
        surfaceVectorFieldCopies_[name] = new surfaceVectorField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, pointScalarField* field)
{
    if (field)
    {
        pointScalarFields_[name] = field;
        pointScalarFieldCopies_[name] = new pointScalarField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, pointVectorField* field)
{
    if (field)
    {
        pointVectorFields_[name] = field;
        pointVectorFieldCopies_[name] = new pointVectorField(*field);
        // TODO: Old time
        // pointVectorFieldsOld_.push_back(const_cast<pointVectorField&>(field->oldTime())));
        // pointVectorFieldCopiesOld_.push_back(new pointVectorField(field->oldTime()));
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, volTensorField* field)
{
    if (field)
    {
        volTensorFields_[name] = field;
        volTensorFieldCopies_[name] = new volTensorField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, surfaceTensorField* field)
{
    if (field)
    {
        surfaceTensorFields_[name] = field;
        surfaceTensorFieldCopies_[name] = new surfaceTensorField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, pointTensorField* field)
{
    if (field)
    {
        pointTensorFields_[name] = field;
        pointTensorFieldCopies_[name] = new pointTensorField(*field);
    }
}

void preciceAdapter::Adapter::addCheckpointField(const std::string& name, volSymmTensorField* field)
{
    if (field)
    {
        volSymmTensorFields_[name] = field;
        volSymmTensorFieldCopies_[name] = new volSymmTensorField(*field);
    }
}


// NOTE: Add here methods to add other object types to checkpoint, if needed.

void preciceAdapter::Adapter::readCheckpoint()
{
    SETUP_TIMER();

    // TODO: To increase efficiency: only the oldTime() fields of the quantities which are used in the time
    //  derivative are necessary. (In general this is only the velocity). Also old information of the mesh
    //  is required.
    //  Therefore, loading the oldTime() and oldTime().oldTime() fields for the other fields can be excluded
    //  for efficiency.
    DEBUG(adapterInfo("Reading a checkpoint..."));

    printFieldCountsDB();

    // Reload the runTime
    reloadCheckpointTime();

    // Reload the meshPoints (if FSI is enabled)
    if (FSIenabled_)
    {
        reloadMeshPoints();
    }


    for (const auto& kv : volScalarFields_)
    {
        const std::string& key = kv.first;
        volScalarField* field = kv.second;
        volScalarField* fieldCopy = volScalarFieldCopies_.at(key);

        // Load the volume field
        *field = *fieldCopy;
        // TODO: Do we need this?
        // field->boundaryField() = fieldCopy->boundaryField();

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    // Reload all the fields of type volVectorField
    for (const auto& kv : volVectorFields_)
    {
        const std::string& key = kv.first;
        volVectorField* field = kv.second;
        volVectorField* fieldCopy = volVectorFieldCopies_.at(key);

        // Load the volume field
        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : surfaceScalarFields_)
    {
        const std::string& key = kv.first;
        Foam::surfaceScalarField* field = kv.second;
        Foam::surfaceScalarField* fieldCopy = surfaceScalarFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : surfaceVectorFields_)
    {
        const std::string& key = kv.first;
        Foam::surfaceVectorField* field = kv.second;
        Foam::surfaceVectorField* fieldCopy = surfaceVectorFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : pointScalarFields_)
    {
        const std::string& key = kv.first;
        Foam::pointScalarField* field = kv.second;
        Foam::pointScalarField* fieldCopy = pointScalarFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : pointVectorFields_)
    {
        const std::string& key = kv.first;
        Foam::pointVectorField* field = kv.second;
        Foam::pointVectorField* fieldCopy = pointVectorFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : volTensorFields_)
    {
        const std::string& key = kv.first;
        Foam::volTensorField* field = kv.second;
        Foam::volTensorField* fieldCopy = volTensorFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : surfaceTensorFields_)
    {
        const std::string& key = kv.first;
        Foam::surfaceTensorField* field = kv.second;
        Foam::surfaceTensorField* fieldCopy = surfaceTensorFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : pointTensorFields_)
    {
        const std::string& key = kv.first;
        Foam::pointTensorField* field = kv.second;
        Foam::pointTensorField* fieldCopy = pointTensorFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    for (const auto& kv : volSymmTensorFields_)
    {
        const std::string& key = kv.first;
        Foam::volSymmTensorField* field = kv.second;
        Foam::volSymmTensorField* fieldCopy = volSymmTensorFieldCopies_.at(key);

        *field = *fieldCopy;

        int nOldTimes = field->nOldTimes();
        if (nOldTimes >= 1)
        {
            field->oldTime() = fieldCopy->oldTime();
        }
        if (nOldTimes == 2)
        {
            field->oldTime().oldTime() = fieldCopy->oldTime().oldTime();
        }
    }

    // NOTE: Add here other field types to read, if needed.

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint was read. Time = " + std::to_string(runTime_.value()));
#endif

    ACCUMULATE_TIMER(timeInCheckpointingRead_);

    return;
}


void preciceAdapter::Adapter::writeCheckpoint()
{
    SETUP_TIMER();

    DEBUG(adapterInfo("Writing a checkpoint..."));

    printFieldCountsDB();

    // Store the runTime
    storeCheckpointTime();

    // Store the meshPoints (if FSI is enabled)
    if (FSIenabled_)
    {
        storeMeshPoints();
    }

    for (const auto& kv : volScalarFields_)
    {
        const std::string& key = kv.first;
        volScalarField* field = kv.second;
        volScalarField* fieldCopy = volScalarFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));
    }

    for (const auto& kv : volVectorFields_)
    {
        const std::string& key = kv.first;
        volVectorField* field = kv.second;
        volVectorField* fieldCopy = volVectorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : volTensorFields_)
    {
        const std::string& key = kv.first;
        volTensorField* field = kv.second;
        volTensorField* fieldCopy = volTensorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : volSymmTensorFields_)
    {
        const std::string& key = kv.first;
        volSymmTensorField* field = kv.second;
        volSymmTensorField* fieldCopy = volSymmTensorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : surfaceScalarFields_)
    {
        const std::string& key = kv.first;
        surfaceScalarField* field = kv.second;
        surfaceScalarField* fieldCopy = surfaceScalarFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : surfaceVectorFields_)
    {
        const std::string& key = kv.first;
        surfaceVectorField* field = kv.second;
        surfaceVectorField* fieldCopy = surfaceVectorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : surfaceTensorFields_)
    {
        const std::string& key = kv.first;
        surfaceTensorField* field = kv.second;
        surfaceTensorField* fieldCopy = surfaceTensorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : pointScalarFields_)
    {
        const std::string& key = kv.first;
        pointScalarField* field = kv.second;
        pointScalarField* fieldCopy = pointScalarFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : pointVectorFields_)
    {
        const std::string& key = kv.first;
        pointVectorField* field = kv.second;
        pointVectorField* fieldCopy = pointVectorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }

    for (const auto& kv : pointTensorFields_)
    {
        const std::string& key = kv.first;
        pointTensorField* field = kv.second;
        pointTensorField* fieldCopy = pointTensorFieldCopies_.at(key);
        *fieldCopy = *field;
        DEBUG(adapterInfo("Checkpointed " + key));

    }
    // NOTE: Add here other types to write, if needed.

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    ACCUMULATE_TIMER(timeInCheckpointingWrite_);

    return;
}

void preciceAdapter::Adapter::readMeshCheckpoint()
{
    DEBUG(adapterInfo("Reading a mesh checkpoint..."));

    // TODO only the meshPhi field is here, which is a surfaceScalarField. The other fields can be removed.
    //  Reload all the fields of type mesh surfaceScalarField
    for (uint i = 0; i < meshSurfaceScalarFields_.size(); i++)
    {
        // Load the volume field
        *(meshSurfaceScalarFields_.at(i)) == *(meshSurfaceScalarFieldCopies_.at(i));

        int nOldTimes(meshSurfaceScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            meshSurfaceScalarFields_.at(i)->oldTime() == meshSurfaceScalarFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            meshSurfaceScalarFields_.at(i)->oldTime().oldTime() == meshSurfaceScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type mesh surfaceVectorField
    for (uint i = 0; i < meshSurfaceVectorFields_.size(); i++)
    {
        // Load the volume field
        *(meshSurfaceVectorFields_.at(i)) == *(meshSurfaceVectorFieldCopies_.at(i));

        int nOldTimes(meshSurfaceVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            meshSurfaceVectorFields_.at(i)->oldTime() == meshSurfaceVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            meshSurfaceVectorFields_.at(i)->oldTime().oldTime() == meshSurfaceVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type mesh volVectorField
    for (uint i = 0; i < meshVolVectorFields_.size(); i++)
    {
        // Load the volume field
        *(meshVolVectorFields_.at(i)) == *(meshVolVectorFieldCopies_.at(i));

        int nOldTimes(meshVolVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            meshVolVectorFields_.at(i)->oldTime() == meshVolVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            meshVolVectorFields_.at(i)->oldTime().oldTime() == meshVolVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh checkpoint was read. Time = " + std::to_string(runTime_.value()));
#endif

    return;
}

void preciceAdapter::Adapter::writeMeshCheckpoint()
{
    DEBUG(adapterInfo("Writing a mesh checkpoint..."));

    // Store all the fields of type mesh surfaceScalar
    for (uint i = 0; i < meshSurfaceScalarFields_.size(); i++)
    {
        *(meshSurfaceScalarFieldCopies_.at(i)) == *(meshSurfaceScalarFields_.at(i));
    }

    // Store all the fields of type mesh surfaceVector
    for (uint i = 0; i < meshSurfaceVectorFields_.size(); i++)
    {
        *(meshSurfaceVectorFieldCopies_.at(i)) == *(meshSurfaceVectorFields_.at(i));
    }

    // Store all the fields of type mesh volVector
    for (uint i = 0; i < meshVolVectorFields_.size(); i++)
    {
        *(meshVolVectorFieldCopies_.at(i)) == *(meshVolVectorFields_.at(i));
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    return;
}

// TODO for the volumes of the mesh, check this part for subcycling.
void preciceAdapter::Adapter::readVolCheckpoint()
{
    DEBUG(adapterInfo("Reading the mesh volumes checkpoint..."));

    // Reload all the fields of type mesh volVectorField::Internal
    for (uint i = 0; i < volScalarInternalFields_.size(); i++)
    {
        // Load the volume field
        *(volScalarInternalFields_.at(i)) = *(volScalarInternalFieldCopies_.at(i));
        // There are no old times for the internal fields.
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh volumes were read. Time = " + std::to_string(runTime_.value()));
#endif

    return;
}

void preciceAdapter::Adapter::writeVolCheckpoint()
{
    DEBUG(adapterInfo("Writing a mesh volumes checkpoint..."));

    // Store all the fields of type mesh volScalarField::Internal
    for (uint i = 0; i < volScalarInternalFields_.size(); i++)
    {
        *(volScalarInternalFieldCopies_.at(i)) = *(volScalarInternalFields_.at(i));
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh volumes checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    return;
}


void preciceAdapter::Adapter::end()
{
    // Throw a warning if the simulation exited before the coupling was complete
    if (NULL != precice_ && isCouplingOngoing())
    {
        adapterInfo("The solver exited before the coupling was complete.", "warning");
    }

    return;
}

void preciceAdapter::Adapter::teardown()
{
    // If the solver interface was not deleted before, delete it now.
    // Normally it should be deleted when isCouplingOngoing() becomes false.
    if (NULL != precice_)
    {
        DEBUG(adapterInfo("Destroying the preCICE solver interface..."));
        delete precice_;
        precice_ = NULL;
    }

    // Delete the preCICE solver interfaces
    if (interfaces_.size() > 0)
    {
        DEBUG(adapterInfo("Deleting the interfaces..."));
        for (uint i = 0; i < interfaces_.size(); i++)
        {
            delete interfaces_.at(i);
        }
        interfaces_.clear();
    }

    // Delete the copied fields for checkpointing
    if (checkpointing_)
    {
        DEBUG(adapterInfo("Deleting the checkpoints... "));

        // Mesh fields
        // meshSurfaceScalar
        for (uint i = 0; i < meshSurfaceScalarFieldCopies_.size(); i++)
        {
            delete meshSurfaceScalarFieldCopies_.at(i);
        }
        meshSurfaceScalarFieldCopies_.clear();

        // meshSurfaceVector
        for (uint i = 0; i < meshSurfaceVectorFieldCopies_.size(); i++)
        {
            delete meshSurfaceVectorFieldCopies_.at(i);
        }
        meshSurfaceVectorFieldCopies_.clear();

        // meshVolVector
        for (uint i = 0; i < meshVolVectorFieldCopies_.size(); i++)
        {
            delete meshVolVectorFieldCopies_.at(i);
        }
        meshVolVectorFieldCopies_.clear();

        // TODO for the internal volume
        //  volScalarInternal
        for (uint i = 0; i < volScalarInternalFieldCopies_.size(); i++)
        {
            delete volScalarInternalFieldCopies_.at(i);
        }
        volScalarInternalFieldCopies_.clear();

        // volScalarField copies
        for (auto& kv : volScalarFieldCopies_)
        {
            delete kv.second; // Delete each dynamic field copy
        }
        volScalarFieldCopies_.clear();

        // volVectorField copies
        for (auto& kv : volVectorFieldCopies_)
        {
            delete kv.second;
        }
        volVectorFieldCopies_.clear();

        // volTensorField copies
        for (auto& kv : volTensorFieldCopies_)
        {
            delete kv.second;
        }
        volTensorFieldCopies_.clear();

        // volSymmTensorField copies
        for (auto& kv : volSymmTensorFieldCopies_)
        {
            delete kv.second;
        }
        volSymmTensorFieldCopies_.clear();

        // surfaceScalarField copies
        for (auto& kv : surfaceScalarFieldCopies_)
        {
            delete kv.second;
        }
        surfaceScalarFieldCopies_.clear();

        // surfaceVectorField copies
        for (auto& kv : surfaceVectorFieldCopies_)
        {
            delete kv.second;
        }
        surfaceVectorFieldCopies_.clear();

        // surfaceTensorField copies
        for (auto& kv : surfaceTensorFieldCopies_)
        {
            delete kv.second;
        }
        surfaceTensorFieldCopies_.clear();

        // pointScalarField copies
        for (auto& kv : pointScalarFieldCopies_)
        {
            delete kv.second;
        }
        pointScalarFieldCopies_.clear();

        // pointVectorField copies
        for (auto& kv : pointVectorFieldCopies_)
        {
            delete kv.second;
        }
        pointVectorFieldCopies_.clear();

        // pointTensorField copies
        for (auto& kv : pointTensorFieldCopies_)
        {
            delete kv.second;
        }
        pointTensorFieldCopies_.clear();


        // NOTE: Add here delete for other types, if needed

        checkpointing_ = false;
    }

    // Delete the CHT module
    if (NULL != CHT_)
    {
        DEBUG(adapterInfo("Destroying the CHT module..."));
        delete CHT_;
        CHT_ = NULL;
    }

    // Delete the FSI module
    if (NULL != FSI_)
    {
        DEBUG(adapterInfo("Destroying the FSI module..."));
        delete FSI_;
        FSI_ = NULL;
    }

    // Delete the FF module
    if (NULL != FF_)
    {
        DEBUG(adapterInfo("Destroying the FF module..."));
        delete FF_;
        FF_ = NULL;
    }

    // NOTE: Delete your new module here

    return;
}

preciceAdapter::Adapter::~Adapter()
{
    teardown();

    TIMING_MODE(
        // Continuing the output started in the destructor of preciceAdapterFunctionObject
        Info << "Time exclusively in the adapter: " << (timeInConfigRead_ + timeInMeshSetup_ + timeInCheckpointingSetup_ + timeInWrite_ + timeInRead_ + timeInCheckpointingWrite_ + timeInCheckpointingRead_).str() << nl;
        Info << "  (S) reading preciceDict:       " << timeInConfigRead_.str() << nl;
        Info << "  (S) constructing preCICE:      " << timeInPreciceConstruct_.str() << nl;
        Info << "  (S) setting up the interfaces: " << timeInMeshSetup_.str() << nl;
        Info << "  (S) setting up checkpointing:  " << timeInCheckpointingSetup_.str() << nl;
        Info << "  (I) writing data:              " << timeInWrite_.str() << nl;
        Info << "  (I) reading data:              " << timeInRead_.str() << nl;
        Info << "  (I) writing checkpoints:       " << timeInCheckpointingWrite_.str() << nl;
        Info << "  (I) reading checkpoints:       " << timeInCheckpointingRead_.str() << nl;
        Info << "  (I) writing OpenFOAM results:  " << timeInWriteResults_.str() << " (at the end of converged time windows)" << nl << nl;
        Info << "Time exclusively in preCICE:     " << (timeInInitialize_ + timeInAdvance_ + timeInFinalize_).str() << nl;
        Info << "  (S) initialize():              " << timeInInitialize_.str() << nl;
        Info << "  (I) advance():                 " << timeInAdvance_.str() << nl;
        Info << "  (I) finalize():                " << timeInFinalize_.str() << nl;
        Info << "  These times include time waiting for other participants." << nl;
        Info << "  See also precice-profiling on the website https://precice.org/tooling-performance-analysis.html." << nl;
        Info << "-------------------------------------------------------------------------------------" << nl;)

    return;
}
