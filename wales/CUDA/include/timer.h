/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/) 
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File timer.h: Timing functionality for CUDA.
 *
 **/

#ifndef TIMER_H
#define TIMER_H

#include <cstring>

#include "error_checking.h"

class Timer
{
	public:
		Timer(const std::string measurementName);
		~Timer();

		void start();
		float stop();
		float elapsed() const;
		void saveMeasurement() const;

		static std::string timerPrefix;

		bool getTimingOn() const;
		void setTimingOn();

	private:
		std::string m_measurementName;

		cudaEvent_t m_start;
		cudaEvent_t m_stop;

		bool m_timingOn;
		bool m_timerRunning;

		float m_accumulatedTime;
};

#endif /* end of include guard: TIMER_H */
