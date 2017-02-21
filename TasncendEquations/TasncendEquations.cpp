// TasncendEquations.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <cmath>
#include <string>
#include <functional>

double Func(double x)
{
	return 8 * cos(x) - x - 6;
}

double DifFunc(double x)
{
	return -8 * sin(x) - 1;
}

double Dichotomy(double Precission, double IntervalStart, double IntervalEnd, std::string FuncStr, std::function<double(double)> Func)
{
	int Dichotomy_max_iter = 100000;
	printf("Dichotomy\n");
	printf("So we trying to find a root of %s in [% 6.4f,% 6.4f]\n", FuncStr.c_str(), IntervalStart, IntervalEnd);
	printf("With precission %e\n", Precission);

	int Steps = 0;

	double ValOnLow, ValOnHigh, ValOnMid;
	double Middle;



	ValOnLow = Func(IntervalStart);
	ValOnHigh = Func(IntervalEnd);

	if (abs(ValOnLow) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalStart);
		return IntervalStart;
	}
	if (abs(ValOnHigh) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalEnd);
		return IntervalEnd;
	}

	if (ValOnLow*ValOnHigh > 0)
	{
		printf("Values in a borders have the same sign. So we ether have no roots in this interval, or an even number of them.\n Anyway, we didn't cut out for this. Terminating\n");
		return 0;
	}


	while (true)
	{
		if (Steps > Dichotomy_max_iter)
		{
			printf("Even after %i iteration, we didn't come up to the root close enough. Terminating\n", Dichotomy_max_iter);
			break;
		}
	
		Middle = (IntervalStart + IntervalEnd) / 2;
		ValOnMid = Func(Middle);

		if (abs(ValOnMid) < Precission)
		{
			printf("So on %ith iteration, we have foud a root %f!\n", Steps, Middle);
			return Middle;
		}

		if (signbit(ValOnMid) == signbit(ValOnLow))
		{//Move the low border
			IntervalStart = Middle;
			ValOnLow = ValOnMid;
		}

		else if (signbit(ValOnMid) == signbit(ValOnHigh))
		{//Move the high border
			IntervalEnd = Middle;
			ValOnHigh = ValOnMid;
		}

		Steps++;
	}

}

double Chordes(double Precission, double IntervalStart, double IntervalEnd, std::string FuncStr, std::function<double(double)> Func)
{
	int Chordes_max_iter = 100000;
	printf("Chordes\n");
	printf("So we trying to find a root of %s in [% 6.4f,% 6.4f]\n", FuncStr.c_str(), IntervalStart, IntervalEnd);
	printf("With precission %e\n", Precission);

	printf("We just assume that everything is fine and input data are viable\n");

	int Steps = 0;

	double ValOnLow, ValOnHigh, ValOnMid;
	double Middle;



	ValOnLow = Func(IntervalStart);
	ValOnHigh = Func(IntervalEnd);

	if (abs(ValOnLow) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalStart);
		return IntervalStart;
	}
	if (abs(ValOnHigh) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalEnd);
		return IntervalEnd;
	}


	if (ValOnLow*ValOnHigh > 0)
	{
		printf("Values in a borders have the same sign. So we ether have no roots in this interval, or an even number of them.\n Anyway, we didn't cut out for this. Terminating\n");
		return 0;
	}

	//if low border is below zero, then high border will be stationary and low border will be moving
	bool LowBorderMoving = signbit(ValOnLow);
	if (LowBorderMoving)
	{
		Middle = IntervalStart;
	}
	else
	{
		Middle = IntervalEnd;
	}

	while (true)
	{
		if (Steps > Chordes_max_iter)
		{
			printf("Even after %i iteration, we didn't come up to the root close enough. Terminating\n", Chordes_max_iter);
			break;
		}

		ValOnMid = Func(Middle);

		if (abs(ValOnMid) < Precission)
		{
			printf("So on %ith iteration, we have foud a root %f!\n", Steps, Middle);
			return Middle;
		}

		if (LowBorderMoving)
		{
			Middle = Middle - ((ValOnMid) / (ValOnMid - ValOnHigh))*(Middle - IntervalEnd);
		}
		else
		{
			Middle = Middle - ((ValOnMid) / (ValOnLow - ValOnMid))*(IntervalStart - Middle);
		}

		Steps++;
	}

}

double Newton(double Precission, double IntervalStart, double IntervalEnd, std::string FuncStr, std::function<double(double)> Func, std::function<double(double)> DifFunc, double Beginning)
{
	int max_iter = 100000;
	printf("Newton\n");
	printf("So we trying to find a root of %s in [% 6.4f,% 6.4f]\n", FuncStr.c_str(), IntervalStart, IntervalEnd);
	printf("With precission %e\n", Precission);

	printf("We just assume that everything is fine and input data are viable\n");
	printf("And we starting from %f\n", Beginning);

	int Steps = 0;

	double ValOnLow, ValOnHigh, ValOnMid;
	double Middle = Beginning;


	ValOnMid = Func(Middle);
	ValOnLow = Func(IntervalStart);
	ValOnHigh = Func(IntervalEnd);

	if ((Middle < IntervalStart) || (Middle > IntervalEnd))
	{
		printf("Starting point is not in the interval\n");
		return 0;
	}

	if (abs(ValOnLow) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalStart);
		return IntervalStart;
	}
	if (abs(ValOnHigh) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalEnd);
		return IntervalEnd;
	}
	

	while (true)
	{
		if (Steps > max_iter)
		{
			printf("Even after %i iteration, we didn't come up to the root close enough. Terminating\n", max_iter);
			break;
		}

		Middle = Middle - Func(Middle) / DifFunc(Middle);
		if (abs(Func(Middle)) < Precission)
		{
			printf("So on %ith iteration, we have foud a root %f!\n", Steps, Middle);
			return Middle;
		}

		if ((Middle < IntervalStart) || (Middle > IntervalEnd))
		{
			printf("Starting point has left the interval. Probably this method is not suit for this function\n");
			return 0;
		}


		Steps++;
	}

}


double NewtonNonDiff(double Precission, double IntervalStart, double IntervalEnd, std::string FuncStr, std::function<double(double)> Func, double Beginning, double PreBeginning)
{
	int max_iter = 100000;
	printf("Newton with no diffs!\n");
	printf("So we trying to find a root of %s in [% 6.4f,% 6.4f]\n", FuncStr.c_str(), IntervalStart, IntervalEnd);
	printf("With precission %e\n", Precission);

	printf("We just assume that everything is fine and input data are viable\n");
	printf("And we starting from %f\n", Beginning);
	printf("With a pre-step in  %f\n", PreBeginning);

	int Steps = 0;

	double ValOnLow, ValOnHigh, ValOnMid;
	double Middle = Beginning;
	double PreMiddle = PreBeginning;

	ValOnMid = Func(Middle);
	ValOnLow = Func(IntervalStart);
	ValOnHigh = Func(IntervalEnd);

	if ((Middle < IntervalStart) || (Middle > IntervalEnd))
	{
		printf("Starting point is not in the interval\n");
		return 0;
	}

	if (abs(ValOnLow) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalStart);
		return IntervalStart;
	}
	if (abs(ValOnHigh) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalEnd);
		return IntervalEnd;
	}

	double TmpLastVal;
	while (true)
	{
		if (Steps > max_iter)
		{
			printf("Even after %i iteration, we didn't come up to the root close enough. Terminating\n", max_iter);
			break;
		}

		TmpLastVal = Middle - Func(Middle) / ((Func(Middle) - Func(PreMiddle)) / (Middle - PreMiddle));
		PreMiddle = Middle;
		Middle = TmpLastVal;



		if (abs(Func(Middle)) < Precission)
		{
			printf("So on %ith iteration, we have foud a root %f!\n", Steps, Middle);
			return Middle;
		}

		if ((Middle < IntervalStart) || (Middle > IntervalEnd))
		{
			printf("Starting point has left the interval. Probably this method is not suit for this function\n");
			return 0;
		}


		Steps++;
	}

}

double Chebishev(double Precission, double IntervalStart, double IntervalEnd, std::string FuncStr, std::function<double(double)> Func, std::function<double(double)> DifFunc, double Beginning, std::function<double(double)> DifDifFunc)
{
	int max_iter = 100000;
	printf("Chebishev\n");
	printf("So we trying to find a root of %s in [% 6.4f,% 6.4f]\n", FuncStr.c_str(), IntervalStart, IntervalEnd);
	printf("With precission %e\n", Precission);

	printf("We just assume that everything is fine and input data are viable\n");
	printf("And we starting from %f\n", Beginning);

	int Steps = 0;

	double ValOnLow, ValOnHigh, ValOnMid;
	double Middle = Beginning;


	ValOnMid = Func(Middle);
	ValOnLow = Func(IntervalStart);
	ValOnHigh = Func(IntervalEnd);

	if ((Middle < IntervalStart) || (Middle > IntervalEnd))
	{
		printf("Starting point is not in the interval\n");
		return 0;
	}

	if (abs(ValOnLow) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalStart);
		return IntervalStart;
	}
	if (abs(ValOnHigh) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalEnd);
		return IntervalEnd;
	}


	while (true)
	{
		if (Steps > max_iter)
		{
			printf("Even after %i iteration, we didn't come up to the root close enough. Terminating\n", max_iter);
			break;
		}

		Middle = Middle - Func(Middle) / DifFunc(Middle) - (DifDifFunc(Middle)*pow(Func(Middle), 2)) / (2 * (pow(DifFunc(Middle), 3)));
		if (abs(Func(Middle)) < Precission)
		{
			printf("So on %ith iteration, we have foud a root %f!\n", Steps, Middle);
			return Middle;
		}

		if ((Middle < IntervalStart) || (Middle > IntervalEnd))
		{
			printf("Starting point has left the interval. Probably this method is not suit for this function\n");
			return 0;
		}


		Steps++;
	}

}


double SimpleIter(double Precission, double IntervalStart, double IntervalEnd, std::string FuncStr, std::function<double(double)> Func, std::string AlterFuncStr, std::function<double(double)> AlterFunc, double Beginning)
{
	int Chordes_max_iter = 100000;
	printf("SimpleIter\n");
	printf("So we trying to find a root of %s in [% 6.4f,% 6.4f]\n", FuncStr.c_str(), IntervalStart, IntervalEnd);
	printf("With precission %e\n", Precission);

	printf("We just assume that everything is fine and input data are viable\n");

	printf("We have %s as iterative func, and we're going to work with that from %f\n", AlterFuncStr.c_str(),Beginning );

	int Steps = 0;

	double ValOnLow, ValOnHigh, ValOnMid;
	double Middle = Beginning;

	if ((Middle < IntervalStart) || (Middle > IntervalEnd))
	{
		printf("Starting point is not in the interval\n");
		return 0;
	}

	ValOnLow = Func(IntervalStart);
	ValOnHigh = Func(IntervalEnd);

	if (abs(ValOnLow) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalStart);
		return IntervalStart;
	}
	if (abs(ValOnHigh) < Precission)
	{
		printf("Surprisingly, border % 6.4f is a root!\n", IntervalEnd);
		return IntervalEnd;
	}



	while (true)
	{
		if (Steps > Chordes_max_iter)
		{
			printf("Even after %i iteration, we didn't come up to the root close enough. Terminating\n", Chordes_max_iter);
			break;
		}

		Middle = AlterFunc(Middle);
		ValOnMid = Func(Middle);



		if (abs(ValOnMid) < Precission)
		{
			printf("So on %ith iteration, we have foud a root %f!\n", Steps, Middle);
			return Middle;
		}

		if ((Middle < IntervalStart) || (Middle > IntervalEnd))
		{
			printf("Starting point has left the interval. Probably this method is not suit for this function\n");
			return 0;
		}

		Steps++;
	}



}


int main()
{
	std::string FunStr = "f(x) = 8 * cos(x) - x - 6";
	std::string FunStrDif = "f(x) = -8 * sin(x) - 1";
	std::function<double(double)> Func = [](double x) {return 8 * cos(x) - x - 6;};
	std::function<double(double)> FuncDif = [](double x) {return -8 * sin(x) - 1;};
	std::function<double(double)> FuncDifDif = [](double x) {return -8 * cos(x);};
	double Preccision = 0.00000000000001;
	double BorderLow = -6;
	double BorderHigh = -4;

	Dichotomy(Preccision, BorderLow, BorderHigh, FunStr, Func);
	Chordes(Preccision, BorderLow, BorderHigh, FunStr, Func);
	Newton(Preccision, BorderLow, BorderHigh, FunStr, Func, FuncDif, -5);
	NewtonNonDiff(Preccision, BorderLow, BorderHigh, FunStr, Func, -4.9, -5);
	Chebishev(Preccision, BorderLow, BorderHigh, FunStr, Func, FuncDif, -5, FuncDifDif);


	std::string AlterFuncStr = "f(x) = x+(3/24)*(8*cos(x)-x-6)";
	std::function<double(double)> AlterFunc = [](double x) {return x + ((double)3/24) * (8 * cos(x) - x - 6);};
	SimpleIter(Preccision, BorderLow, BorderHigh, FunStr, Func, AlterFuncStr, AlterFunc, -5 );
    return 0;
}

