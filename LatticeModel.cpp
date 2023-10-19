#include "LatticeModel.h"

LatticeModel::LatticeModel()
	: VelocitySet(
		{ 8.0 / 27.0,      //0
		  2.0 / 27.0,      //1 
		  2.0 / 27.0,      //2
		  2.0 / 27.0,      //3
		  2.0 / 27.0,      //4
		  2.0 / 27.0,      //5
		  2.0 / 27.0,      //6
		  1.0 / 54.0,      //7
		  1.0 / 54.0,      //8
		  1.0 / 54.0,      //9
		  1.0 / 54.0,      //10
		  1.0 / 54.0,      //11
		  1.0 / 54.0,      //12
		  1.0 / 54.0,      //13
		  1.0 / 54.0,      //14
		  1.0 / 54.0,      //15
		  1.0 / 54.0,      //16
		  1.0 / 54.0,      //17
		  1.0 / 54.0,      //18
		  1.0 / 216.0,     //19
		  1.0 / 216.0,     //20
		  1.0 / 216.0,     //21
		  1.0 / 216.0,     //22
		  1.0 / 216.0,     //23
		  1.0 / 216.0,     //24
		  1.0 / 216.0,     //25
		  1.0 / 216.0,     //26

        },
		{ { 0,  0, 0},     //0
		  { 1,  0, 0},     //1
		  {-1,  0, 0},     //2
		  { 0, -1, 0},     //3
		  { 0,  1, 0},     //4 
		  { 0,  0, 1},     //5 
		  { 0,  0,-1},     //6
		  { 1, -1, 0},     //7   
		  {-1, -1, 0},     //8
		  { 1,  1, 0},     //9
		  {-1,  1, 0},     //10
		  { 1,  0, 1},     //11
		  {-1,  0, 1},     //12
		  { 1,  0,-1},     //13
		  {-1,  0,-1},     //14
		  { 0, -1, 1},     //15
		  { 0,  1, 1},     //16
		  { 0, -1,-1},     //17
		  { 0,  1,-1},     //18
		  { 1, -1, 1},     //19
		  {-1, -1, 1},     //20
		  { 1,  1, 1},     //21
		  {-1,  1, 1},     //22
		  { 1, -1,-1},     //23
		  {-1, -1,-1},     //24
		  { 1,  1,-1},     //25
		  {-1,  1,-1},     //26


		},
		1.0 / 3.0
	)
{
}


LatticeModel::~LatticeModel()
{
}

int LatticeModel::OppositeDirection(int dir) const
{
	switch (dir)
	{
	case 0:
		return 0;
	case 1:
		return 2;
	case 2:
		return 1;
	case 3:
		return 4;
	case 4:
		return 3;
	case 5:
		return 6;
	case 6:
		return 5;
	case 7:
		return 10;
	case 8:
		return 9;
	case 9:
		return 8;
	case 10:
		return 7;
	case 11:
		return 14;
	case 12:
		return 13;
	case 13:
		return 12;
	case 14:
		return 11;
	case 15:
		return 18;
	case 16:
		return 17;
	case 17:
		return 16;
	case 18:
		return 15;
	case 19:
		return 26;
	case 20:
		return 25;
	case 21:
		return 24;
	case 22:
		return 23;
	case 23:
		return 22;
	case 24:
		return 21;
	case 25:
		return 20;
	case 26:
		return 19;

	}
}