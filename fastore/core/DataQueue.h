/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DATAQUEUE
#define H_DATAQUEUE

#include "Globals.h"

#include <queue>

#include "Thread.h"


/**
 * Templatied data queue used in multi-threaded processing
 *
 */
template <class _TDataType>
class TDataQueue
{
	typedef _TDataType DataType;
	typedef std::pair<int64, DataType*> part_pair;
	typedef std::deque<part_pair> part_queue;

public:
	static const uint32 DefaultMaxPartNum = 64;

	TDataQueue(uint32 maxPartNum_ = DefaultMaxPartNum, uint32 threadNum_ = 1)
		:	threadNum(threadNum_)
		,	maxPartNum(maxPartNum_)
		,	currentThreadMask(0)
	{
		ASSERT(maxPartNum_ > 0);
		ASSERT(threadNum_ >= 1);
		ASSERT(threadNum_ < 64);

		completedThreadMask = ((uint64)1 << threadNum) - 1;
	}

	~TDataQueue()
	{}

	bool IsEmpty()
	{
		return parts.empty();
	}

	bool IsCompleted()
	{
		return parts.empty() && currentThreadMask == completedThreadMask;
	}

	void SetCompleted()
	{
		mt::lock_guard<mt::mutex> lock(mutex);

		ASSERT(currentThreadMask != completedThreadMask);
		currentThreadMask = (currentThreadMask << 1) | 1;

		queueEmptyCondition.notify_all();
	}

	void Push(int64 partId_, const DataType* part_)
	{
		mt::unique_lock<mt::mutex> lock(mutex);

		while (parts.size() > maxPartNum && partId_ > parts.front().first)
			queueFullCondition.wait(lock);

		parts.push_back(std::make_pair(partId_, (DataType*)part_));
		if(parts.size() > 1)
		{
			std::sort(parts.begin(), parts.end(), sortByPartId);
		}

		queueEmptyCondition.notify_one();
	}

	bool Pop(int64 &partId_, DataType* &part_)
	{
		mt::unique_lock<mt::mutex> lock(mutex);

		while ((parts.size() == 0) && currentThreadMask != completedThreadMask)
			queueEmptyCondition.wait(lock);

		if (parts.size() != 0)
		{
			partId_ = parts.front().first;
			part_ = parts.front().second;
			parts.pop_front();
			queueFullCondition.notify_one();
			return true;
		}

		// assure this is impossible
		ASSERT(currentThreadMask == completedThreadMask);
		ASSERT(parts.size() == 0);
		return false;
	}

	void Reset()
	{
		ASSERT(currentThreadMask == completedThreadMask);
		ASSERT(parts.size() == 0);
		currentThreadMask = 0;
	}

private:
	const uint32 threadNum;
	const uint32 maxPartNum;
	uint64 completedThreadMask;
	uint64 currentThreadMask;
	part_queue parts;

	mt::mutex mutex;
	mt::condition_variable queueFullCondition;
	mt::condition_variable queueEmptyCondition;

	inline static bool sortByPartId(const part_pair& p1_, const part_pair& p2_)
	{
		return p1_.first < p2_.first;
	}

};

#endif // H_DATAQUEUE
