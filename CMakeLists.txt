cmake_minimum_required(VERSION 3.5)

project(CEOS)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f03 -free -implicitnone -mkl -fopenmp")
    set(bounds "-check bounds")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

add_executable(
    ceos
        SourceCode/Tools/ModIO.f90
        SourceCode/Tools/SubError.f90
        SourceCode/Tools/ModTimer.f90
        SourceCode/Tools/ModCharacter.f90
        SourceCode/Tools/ModStatus.f90
        SourceCode/Tools/ModTools.f90
        SourceCode/Tools/MathTools/SparseMatrixRoutines/ModSparseVectorRoutines.f90
        SourceCode/Tools/MathTools/SparseMatrixRoutines/ModSparseMatrixRoutines.f90
        SourceCode/Tools/MathTools/SparseMatrixRoutines/ModGlobalSparseMatrix.f90
        SourceCode/Tools/MathTools/ModTensorAlgebra.f90
        SourceCode/Tools/MathTools/ModMathRoutines.f90
        SourceCode/Tools/MathTools/ModVoigtNotation.f90
        SourceCode/Tools/MathTools/ModMathRoutines_NEW.f90
        SourceCode/Tools/ModParser.f90
        SourceCode/Tools/MathTools/LinearSolver/ModLinearSolver.f90
        SourceCode/Tools/MathTools/LinearSolver/ModPardisoSolver.f90
        SourceCode/Tools/MathTools/LinearSolver/ModFullLinearSolver.f90
        SourceCode/Tools/MathTools/LinearSolver/ModLinearSolverLibrary.f90
        SourceCode/Tools/MathTools/ModNonLinearSystemOfEquations.f90
        SourceCode/Tools/MathTools/NonlinearSolver/ModNonLinearSolver.f90
        SourceCode/Tools/MathTools/NonlinearSolver/ModNewtonRaphsonFull.f90
        SourceCode/Tools/MathTools/NonlinearSolver/ModNonLinearSolverLibrary.f90
        SourceCode/Tools/ModContinuumMechanics.f90
        SourceCode/FEM/Modules/ModNodes.f90
        SourceCode/FEM/Modules/ModAnalysis.f90
        SourceCode/ConstitutiveLibrary/ModConstitutiveModel.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModGeneralizedHookesLaw.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModJ2Plasticity.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModNeoHookeanIsochoric.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModViscoelasticFiber.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModViscoelasticMatrix.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModStVenantKirchhoff.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModViscoelasticMatrixFiber.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModNeoHookean.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModNeoHookeanQ1P0.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModHyperelasticQ1P0.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModHyperelasticTransIsoComp.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModHyperelasticTransIso.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModCompressibleNeoHookean.f90
        SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModConstitutiveModelLibrary.f90
        SourceCode/FEM/Elements/ElementMethods/ModElement.f90
        SourceCode/FEM/Elements/ElementLibrary/ModElementQuad4.f90
        SourceCode/FEM/Elements/ElementLibrary/ModElementTri3.f90
        SourceCode/FEM/Elements/ElementLibrary/ModElementTetra10.f90
        SourceCode/FEM/Elements/ElementLibrary/ModElementTetra4.f90
        SourceCode/FEM/Elements/ElementLibrary/ModElementHexa8.f90
        SourceCode/FEM/Elements/ElementLibrary/ModElementLibrary.f90
        SourceCode/FEM/Elements/ElementMethods/SubElementConstructor.f90
        SourceCode/ConstitutiveLibrary/SubMaterialConstructor.f90
        SourceCode/FEM/Modules/ModuleInterfaces.f90
        SourceCode/FEM/BoundaryConditions/ModLoadHistoryData.f90
        SourceCode/FEM/BoundaryConditions/ModBoundaryConditions.f90
        SourceCode/FEM/Modules/ModFEMSystemOfEquations.f90
        SourceCode/FEM/Subroutines/SubAssembleGlobalMatrixUpperTriangular.f90
        SourceCode/FEM/Subroutines/SubUpdateMeshCoordinates.f90
        SourceCode/FEM/Subroutines/SubArgumentHandler.f90
        SourceCode/FEM/Subroutines/SubAssembleGlobalMatrix.f90
        SourceCode/FEM/Subroutines/SubTangentStiffnessMatrix.f90
        SourceCode/FEM/Subroutines/SubInternalForce.f90
        SourceCode/FEM/Subroutines/SubSolveConstitutiveModel.f90
        SourceCode/FEM/Subroutines/SubExternalForce.f90
        SourceCode/FEM/Subroutines/SubPrecribedDisplacement.f90
        SourceCode/MultiScale/ModMultiscaleBoundaryConditions.f90
        SourceCode/FEM/Input/ModReadInputFile.f90
        SourceCode/FEM/ModFEMAnalysis.f90
        SourceCode/MultiScale/ModMultiscaleAnalysis.f90
        SourceCode/ModAnalysisManager.f90
        SourceCode/FEM/Input/SubAnalyzeLoadCaseTables.f90
        SourceCode/FEM/PostProcessing/ModProbe.f90
        SourceCode/FEM/PostProcessing/ModPostProcessors.f90
        SourceCode/FEM/PostProcessing/ModGiD.f90
        SourceCode/FEM/PostProcessing/ModHyperView.f90
        SourceCode/FEM/PostProcessing/ModExportResultFile.f90
        SourceCode/MAIN.f90
)
