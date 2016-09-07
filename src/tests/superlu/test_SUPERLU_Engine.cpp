#include "slu_ddefs.h"
#include <vector>
#include "chrono_superlu/ChSuperLUEngine.h"
#include <core/ChCSR3Matrix.h>

using namespace chrono;


int test_pure_vectors_input()
{
	std::cout << "/******* Test simple solution using STL vectors *******/" << std::endl;
	std::vector<double> values = { 19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18 };
	std::vector<int> rowIndex = { 0,1,4,1,2,4,0,2,0,3,3,4 };
	std::vector<int> colIndex = { 0,3,6,8,10,12 };

	int rows = colIndex.size() - 1;

	std::vector<double> rhs(rows, 1.0);
	std::vector<double> sol(rows, 0.0);


	ChSuperLUEngine m_engine;

	m_engine.SetMatrix(rows, values.data(), rowIndex.data(), colIndex.data());
	m_engine.SetSolutionVector(sol.data());
	m_engine.SetRhsVector(rhs.data());

	m_engine.SuperLUCall(13);

	for (auto row_sel = 0; row_sel<sol.size(); ++row_sel)
	{
		std::cout << sol[row_sel] << " ";
	}
	std::cout << std::endl << std::endl << std::endl;

	return 0;
}

int test_CSR_input()
{
	std::cout << "/******* Test simple solution using CSR matrix *******/" << std::endl;
	int rows = 5;
	ChCSR3Matrix matCSR3(rows, rows, false);

	matCSR3.SetElement(0, 0, 19);
	matCSR3.SetElement(1, 0, 12);
	matCSR3.SetElement(4, 0, 12);

	matCSR3.SetElement(1, 1, 21);
	matCSR3.SetElement(2, 1, 12);
	matCSR3.SetElement(4, 1, 12);

	matCSR3.SetElement(0, 2, 21);
	matCSR3.SetElement(2, 2, 16);

	matCSR3.SetElement(0, 3, 21);
	matCSR3.SetElement(3, 3, 5);

	matCSR3.SetElement(3, 4, 21);
	matCSR3.SetElement(4, 4, 18);

	matCSR3.Compress();


	GetLog() << "Input matrix\n";
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < rows; j++)
			printf("%f ", matCSR3.GetElement(i, j));

		printf("\n");
	}

	

	std::vector<double> rhs(rows, 1.0);
	std::vector<double> sol(rows, 0.0);


	ChSuperLUEngine m_engine;

	m_engine.SetMatrix(matCSR3);
	m_engine.SetSolutionVector(sol.data());
	m_engine.SetRhsVector(rhs.data());

	m_engine.SuperLUCall(13, true);

	for (auto row_sel = 0; row_sel<sol.size(); ++row_sel)
		std::cout << sol[row_sel] << " ";
	std::cout << std::endl << std::endl << std::endl;

	return 0;
}

int test_SetupSolve()
{
	std::cout << " /******* Test Factorization + Solution phase *******/" << std::endl;
	// initialize matrix
	int rows = 5;
	ChCSR3Matrix matCSR3(rows, rows, false);

	matCSR3.SetElement(0, 0, 19);
	matCSR3.SetElement(1, 0, 12);
	matCSR3.SetElement(4, 0, 12);

	matCSR3.SetElement(1, 1, 21);
	matCSR3.SetElement(2, 1, 12);
	matCSR3.SetElement(4, 1, 12);

	matCSR3.SetElement(0, 2, 21);
	matCSR3.SetElement(2, 2, 16);

	matCSR3.SetElement(0, 3, 21);
	matCSR3.SetElement(3, 3, 5);

	matCSR3.SetElement(3, 4, 21);
	matCSR3.SetElement(4, 4, 18);

	matCSR3.Compress();

	std::vector<double> rhs(rows, 1.0);
	std::vector<double> sol(rows, 0.0);


	ChSuperLUEngine m_engine;
	// call for factorization
	m_engine.SetMatrix(matCSR3);
	m_engine.SuperLUCall(12, true);

	// call for first solve
	m_engine.SetRhsVector(rhs.data());
	m_engine.SetSolutionVector(sol.data());
	m_engine.SuperLUCall(33, true);

	for (auto row_sel = 0; row_sel<sol.size(); ++row_sel)
		std::cout << sol[row_sel] << " ";
	std::cout << std::endl << std::endl;

	// call for second solve
	std::fill(rhs.begin(), rhs.end(), 2.0);
	m_engine.SuperLUCall(33, true);

	for (auto row_sel = 0; row_sel<sol.size(); ++row_sel)
		std::cout << sol[row_sel] << " ";
	std::cout << std::endl << std::endl << std::endl;
	return 0;
}

int main(int argc, char *argv[])
{

	int test1 = test_pure_vectors_input();
	int test2 = test_CSR_input();
	int test3 = test_SetupSolve();

	

	


	getchar();
	
	return test1 || test2 || test3;
}
