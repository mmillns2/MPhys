#pragma once


#include<chrono>


// To use, instantiate a Timer at start of code and call Timer::elapsed() at end of code.


class Timer
{
public:
	void reset() { m_beg = Clock::now(); }
	double elapsed() const { return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count(); }

private:
	using Clock = std::chrono::steady_clock;
	using Second = std::chrono::duration<double, std::ratio<1>>;
 	std::chrono::time_point<Clock> m_beg { Clock::now() };
};
