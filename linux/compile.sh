ifort -mkl -fopenmp -c ../SourceCode/Tools/ModIO.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/SubError.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/ModTimer.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/ModCharacter.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/ModStatus.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/ModTools.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/SparseMatrixRoutines/ModSparseVectorRoutines.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/SparseMatrixRoutines/ModSparseMatrixRoutines.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/SparseMatrixRoutines/ModGlobalSparseMatrix.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/ModTensorAlgebra.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/ModMathRoutines.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/ModVoigtNotation.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/ModMathRoutines_NEW.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/ModParser.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/LinearSolver/ModLinearSolver.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/LinearSolver/ModPardisoSolver.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/LinearSolver/ModFullLinearSolver.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/LinearSolver/ModLinearSolverLibrary.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/ModNonLinearSystemOfEquations.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/NonlinearSolver/ModNonLinearSolver.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/NonlinearSolver/ModNewtonRaphsonFull.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/MathTools/NonlinearSolver/ModNonLinearSolverLibrary.f90
ifort -mkl -fopenmp -c ../SourceCode/Tools/ModContinuumMechanics.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Modules/ModNodes.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Modules/ModAnalysis.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ModConstitutiveModel.f90
ifort -mkl -fopenmp -c "../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModGeneralizedHookesLaw .f90"
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModJ2Plasticity.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModNeoHookeanIsochoric.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModViscoelasticFiber.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModViscoelasticMatrix.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModStVenantKirchhoff.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModViscoelasticMatrixFiber.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModNeoHookean.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModNeoHookeanQ1P0.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModHyperelasticQ1P0.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModHyperelasticTransIsoComp.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModHyperelasticTransIso.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModCompressibleNeoHookean.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/ConstitutiveModelsLibrary/ModConstitutiveModelLibrary.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementMethods/ModElement.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementLibrary/ModElementQuad4.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementLibrary/ModElementTri3.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementLibrary/ModElementTetra10.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementLibrary/ModElementTetra4.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementLibrary/ModElementHexa8.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementLibrary/ModElementLibrary.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Elements/ElementMethods/SubElementConstructor.f90
ifort -mkl -fopenmp -c ../SourceCode/ConstitutiveLibrary/SubMaterialConstructor.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Modules/ModuleInterfaces.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/BoundaryConditions/ModLoadHistoryData.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/BoundaryConditions/ModBoundaryConditions.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Modules/ModFEMSystemOfEquations.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubAssembleGlobalMatrixUpperTriangular.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubUpdateMeshCoordinates.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubArgumentHandler.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubAssembleGlobalMatrix.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubTangentStiffnessMatrix.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubInternalForce.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubSolveConstitutiveModel.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubExternalForce.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Subroutines/SubPrecribedDisplacement.f90
ifort -mkl -fopenmp -c ../SourceCode/MultiScale/ModMultiscaleBoundaryConditions.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Input/ModReadInputFile.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/ModFEMAnalysis.f90
ifort -mkl -fopenmp -c ../SourceCode/MultiScale/ModMultiscaleAnalysis.f90
ifort -mkl -fopenmp -c ../SourceCode/ModAnalysisManager.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/Input/SubAnalyzeLoadCaseTables.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/PostProcessing/ModProbe.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/PostProcessing/ModPostProcessors.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/PostProcessing/ModGiD.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/PostProcessing/ModHyperView.f90
ifort -mkl -fopenmp -c ../SourceCode/FEM/PostProcessing/ModExportResultFile.f90

ifort -mkl *.o ../SourceCode/MAIN.f90 -o ceos

echo "Done"
