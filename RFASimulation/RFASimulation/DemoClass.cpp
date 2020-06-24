
#include "DemoClass.h";

using namespace RFASimulation;


DemoClass::DemoClass(int myInt)
{
    DemoClass::m_myInt = myInt;
}


DemoClass::~DemoClass(void) {};

int DemoClass::GetMyInt()
{
    return m_myInt;
}
