#include "slu_ddefs.h"
#include <vector>
#include "chrono_superlu/ChSuperLUEngine.h"

using namespace chrono;

int main(int argc, char *argv[])
{

	std::vector<double> values = { 19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18};
	std::vector<int> rowIndex = { 0,1,4,1,2,4,0,2,0,3,3,4 };
	std::vector<int> colIndex = { 0,3,6,8,10,12 };

	int rows = colIndex.size()-1;

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

	getchar();


}
