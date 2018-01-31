/**
 *
 * File printing.cu: Implementation of class Printing.
 *
 **/

#include "printing.h"

	Printing::Printing(const std::string fileName)
: m_printingOn(false)
{
	m_fileHandle.open(fileName.c_str(), std::ios_base::app);
}

Printing::~Printing()
{
	m_fileHandle.close();
}

bool Printing::getPrintingOn() const
{
	return m_printingOn;
}

void Printing::setPrintingOn()
{
	m_printingOn = true;
}

std::ofstream &Printing::getFileHandle()
{
	return m_fileHandle;
}
