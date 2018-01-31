/**
 *
 * File printing.h: Printing class for outputting debug information. 
 *
 **/

#ifndef PRINTING_H
#define PRINTING_H

#include <cstring>
#include <fstream>

class Printing
{
	public:
		Printing(const std::string fileName);
		~Printing();

		bool getPrintingOn() const;
		void setPrintingOn();

		std::ofstream &getFileHandle();

	private:
		std::ofstream m_fileHandle;

		bool m_printingOn; // True if debug printing is turned on. 
};

#endif /* end of include guard: PRINTING_H */
