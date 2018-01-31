/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/) 
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File timer_cuda.cu: Implementation of class Timer.
 *
 **/

#include <fstream>

#include "timer.h"

	Timer::Timer(const std::string measurementName)
	: m_measurementName(measurementName)
	, m_timerRunning(false)
	, m_accumulatedTime(0.0f)
	  , m_timingOn(false)
{
	CudaSafeCall( cudaEventCreate(&m_start) );
	CudaSafeCall( cudaEventCreate(&m_stop)  );
}

Timer::~Timer()
{
	if (m_timerRunning) {
		stop();
		saveMeasurement();
	}

	CudaSafeCall( cudaEventDestroy(m_start) );
	CudaSafeCall( cudaEventDestroy(m_stop)  );
}

bool Timer::getTimingOn() const
{
	return m_timingOn;
}

void Timer::setTimingOn()
{
	m_timingOn = true;
}

void Timer::start()
{
	m_timerRunning = true;
	CudaSafeCall( cudaEventRecord(m_start, 0) );
}

float Timer::stop()
{
	CudaSafeCall( cudaEventRecord(m_stop, 0) );
	CudaSafeCall( cudaEventSynchronize(m_stop) );
	m_timerRunning = false;

	float elapsedTime;
	CudaSafeCall( cudaEventElapsedTime(&elapsedTime, m_start, m_stop) );

	m_accumulatedTime += elapsedTime;

	return elapsed();
}

void Timer::saveMeasurement() const
{
	std::string filename(timerPrefix);
	filename.append(m_measurementName);
	filename.append(".txt");

	std::ofstream stream;
	stream.open(filename.c_str(), std::ios_base::app);
	stream << elapsed() << std::endl;
	stream.close();
}

float Timer::elapsed() const
{
	return m_accumulatedTime;
}

std::string Timer::timerPrefix = "";

