/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"
#include "FastqParser.h"

#include "../qvz/Stats.h"	// for updating stats definition


bool SingleDnaRecordParser::ReadLine(uchar *str_, uint32& len_, uint32& size_)
{
	uint32 i = 0;
	for (; ;)
	{
		int32 c = Getc();
		if (c == -1)
			break;

		if (c != '\n' && c != '\r')
		{
			if (i >= size_)
			{
				extend_string(str_, size_);
			}
			str_[i++] = (uchar)c;
		}
		else
		{
			if (c == '\r' && Peekc() == '\n')	// case of CR LF
				Skipc();

			if (i > 0)
				break;
		}
	}
	str_[i] = 0;
	len_ = i;
	return i > 0;
}


uint32 SingleDnaRecordParser::SkipLine()
{
	uint32 len = 0;
	for (;;)
	{
		int32 c = Getc();
		if (c == -1)
			break;

		if (c != '\n' && c != '\r')
		{
			len++;
		}
		else
		{
			if (c == '\r' && Peekc() == '\n')	// case of CR LF
				Skipc();

			break;
		}
	}
	return len;
}


bool SingleDnaRecordParser::ReadNextRecord(FastqRecord& rec_)
{
	if (memoryPos >= memorySize)
		return false;

	char* title = (char*)(memory + memoryPos);
	uint32 titleLen = SkipLine();
	if (titleLen == 0 || title[0] != '@')
		return false;

	char* seq = (char*)(memory + memoryPos);
	uint32 seqLen = SkipLine();
	ASSERT(seqLen < FastqRecord::MaxSeqLen);

	uint16 plen = SkipLine();
	if (plen == 0)
		return false;

	//char* qua = (char*)(memory + memoryPos);
	uint16 qlen = SkipLine();
	if (qlen != seqLen)
		return false;

	rec_.seq = seq;
	rec_.seqLen = seqLen;

	stats->Update(rec_);

	return true;
}


void SingleDnaRecordParser::WriteNextRecord(const FastqRecord& rec_)
{
	ASSERT(rec_.seqLen > 0);
	if (memoryPos + rec_.seqLen + 1 > memorySize)
		ExtendBuffer(memoryPos + rec_.seqLen + 1);

	std::copy(rec_.seq, rec_.seq + rec_.seqLen, memory + memoryPos);
	memoryPos += rec_.seqLen;

	memory[memoryPos++] = '\n';
}




bool SingleFastqRecordParser::ReadNextRecord(FastqRecord& rec_)
{
	if (memoryPos >= memorySize)
		return false;

	char* title = (char*)(memory + memoryPos);
	uint32 titleLen = SkipLine();
	if (titleLen == 0 || title[0] != '@')
		return false;

	char* seq = (char*)(memory + memoryPos);
	uint32 seqLen = SkipLine();
	ASSERT(seqLen < FastqRecord::MaxSeqLen);

	uint16 plen = SkipLine();
	if (plen == 0)
		return false;

	char* qua = (char*)(memory + memoryPos);
	uint16 qlen = SkipLine();
	if (qlen != seqLen)
		return false;

	rec_.seq = seq;
	rec_.seqLen = seqLen;
	rec_.qua = qua;
	if (useHeaders)
	{
		rec_.head = title;

		if (!keepComments)
		{
			for (rec_.headLen = 0 ; rec_.headLen < titleLen; rec_.headLen++)
			{
				if (rec_.head[rec_.headLen] == ' ')
					break;
			}
		}
		else
		{
			rec_.headLen = titleLen;
		}
	}

	stats->Update(rec_);

	return true;
}


void SingleFastqRecordParser::WriteNextRecord(const FastqRecord& rec_)
{
	ASSERT(rec_.seqLen > 0);
	const uint64 approxSize = (rec_.seqLen * 2) + rec_.headLen + 5;		// 5: 4 newline + '+'
	if (memoryPos + approxSize > memorySize)
		ExtendBuffer((((memoryPos + approxSize) / 4096) + 1) * 4096);


	// store header
	//
	/*
	if (readsHaveHeaders)
	{
		ASSERT(rec_.head != NULL);
		ASSERT(rec_.headLen > 0);

		std::copy(rec_.head, rec_.head + rec_.headLen, memory + memoryPos);

		memoryPos += rec_.headLen;
	}
	else
	{
		uchar ibuf[32];
		uint32 headLen = to_string(ibuf, recIdx++);		// TODO: this can be easily optimized

		std::copy(headerPrefix.begin(), headerPrefix.end(), memory + memoryPos);
		std::copy(ibuf, ibuf + headLen, memory + memoryPos + headerPrefix.size());

		memoryPos += headerPrefix.size() + headLen;
	}
	*/

	ASSERT(rec_.head != NULL);
	ASSERT(rec_.headLen > 0);
	std::copy(rec_.head, rec_.head + rec_.headLen, memory + memoryPos);
	memoryPos += rec_.headLen;

	memory[memoryPos++] = '\n';


	// store sequence
	//
	std::copy(rec_.seq, rec_.seq + rec_.seqLen, memory + memoryPos);
	memoryPos += rec_.seqLen;
	memory[memoryPos++] = '\n';


	// store comment
	//
	memory[memoryPos++] = '+';
	memory[memoryPos++] = '\n';


	// store quality
	//
	std::copy(rec_.qua, rec_.qua + rec_.seqLen, memory + memoryPos);
	memoryPos += rec_.seqLen;
	memory[memoryPos++] = '\n';
}


void SingleDnaRecordParser::StartParsing(DataChunk &chunk_, ParserMode mode_, FastqRawBlockStats* stats_)
{
	buf = &chunk_.data;
	memory = buf->Pointer();
	memoryPos = 0;

	if (mode_ == ParseRead)
	{
		ASSERT(stats_ != NULL);
		memorySize = chunk_.size;
		stats = stats_;
	}
	else
	{
		memorySize = chunk_.data.Size();
	}
}

uint64 SingleDnaRecordParser::FinishParsing(ParserMode mode_)
{
	if (mode_ == ParseRead)
		return memorySize - skippedBytes;
	else
		return memoryPos;
}


// full FASTQ parsers
//
uint64 FastqRecordsParserSE::ParseTo(const std::vector<FastqRecord> &reads_, IFastqChunkCollection &chunk_, uint64 recStartIdx_)
{
	ASSERT(!chunk_.chunks.empty());

	SingleFastqRecordParser parser(useHeaders);
	parser.StartParsing(*chunk_.chunks[0], SingleDnaRecordParser::ParseWrite);

	FastqRecordBuffer rcBuf;

	std::array<char, FastqRecord::MaxTagLen> headBuf;
	if (!useHeaders)
		std::copy(autoHeaderPrefix.begin(), autoHeaderPrefix.end(), headBuf.begin());

	for (const auto& rec : reads_)
	{
		if (rec.IsReadReverse())
		{
			rec.ComputeRC(rcBuf);
			if (!useHeaders)
			{
				rcBuf.head = headBuf.data();
				rcBuf.headLen = autoHeaderPrefix.size();
				rcBuf.headLen += to_string(headBuf.begin() + autoHeaderPrefix.size(), recStartIdx_);		// <library>.<num>
			}
			else
			{
				rcBuf.head = rec.head;
				rcBuf.headLen = rec.headLen;
			}

			parser.WriteNextRecord(rcBuf);
		}
		else
		{
			if (!useHeaders)
			{
				FastqRecord r = rec;

				r.head = headBuf.data();
				r.headLen = autoHeaderPrefix.size();
				r.headLen += to_string(headBuf.begin() + autoHeaderPrefix.size(), recStartIdx_);		// <library>.<num>

				parser.WriteNextRecord(r);
			}
			else
			{
				parser.WriteNextRecord(rec);
			}
		}
		recStartIdx_++;
	}

	chunk_.chunks[0]->size = parser.FinishParsing(SingleDnaRecordParser::ParseWrite);
	return chunk_.chunks[0]->size;
}


uint64 FastqRecordsParserSE::ParseFrom(IFastqChunkCollection& chunk_, std::vector<FastqRecord>& records_, FastqRawBlockStats& stats_, bool keepComments_)
{
	ASSERT(!chunk_.chunks.empty());

	records_.clear();
	stats_.Clear();

#if EXTRA_MEM_OPT
	records_.shrink_to_fit();
#endif

	// we're losing here const qualifier -- we need to redesign underlying parser
	SingleFastqRecordParser parser(useHeaders, keepComments_);
	parser.StartParsing(*chunk_.chunks[0], SingleDnaRecordParser::ParseRead, &stats_);


	// INFO: while parsing we can also move and copy data inside buffer for better locality
	//
	FastqRecord rec;
	while (parser.ReadNextRecord(rec))
	{
		rec.SetReadReverse(false);
		rec.minimPos = 0;
		records_.emplace_back(rec);
	}
	ASSERT(records_.size() > 0);

	return parser.FinishParsing(SingleDnaRecordParser::ParseRead);
}


uint64 FastqRecordsParserPE::ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_)
{
	ASSERT(chunk_.chunks.size() >= 2);

	SingleFastqRecordParser parser_1(useHeaders), parser_2(useHeaders);

	parser_1.StartParsing(*chunk_.chunks[0], SingleDnaRecordParser::ParseWrite);
	parser_2.StartParsing(*chunk_.chunks[1], SingleDnaRecordParser::ParseWrite);

	FastqRecordBuffer recBuf;

	std::array<char, FastqRecord::MaxTagLen> headBuf;

	if (!useHeaders)
		std::copy(autoHeaderPrefix.begin(), autoHeaderPrefix.end(), headBuf.begin());

	for (const FastqRecord& rec : reads_)
	{
		if (useHeaders)
		{
			// here we set the proper headers
			//
			std::copy(rec.head, rec.head + rec.headLen, headBuf.begin());


			// check whether the /1 /2 fields are present
			//
			if (peFieldIdx != 0)
			{
				// find position of the /2 token
				const std::string separators = FastqRawBlockStats::HeaderStats::Separators();
				ASSERT(peFieldIdx != 0);
				uint32 pairTokenPos = 0;

				uint32 fidx = 0;
				for (uint32 i = 0; i <= rec.headLen; ++i)
				{
					if (!std::count(separators.begin(), separators.end(), rec.head[i]) && (i != rec.headLen))
						continue;

					fidx++;

					if (fidx == peFieldIdx)
					{
						pairTokenPos = i+1;
						break;
					}
				}

				ASSERT(rec.head[pairTokenPos] == '1');
				headBuf[pairTokenPos] = '2';
			}


			if (rec.IsReadReverse())
			{
				rec.ComputeRC(recBuf);
				FastqRecord rAux = recBuf.GetPair();

				if (rec.IsPairSwapped())
				{
					rAux.head = rec.head;
					rAux.headLen = rec.headLen;
					parser_1.WriteNextRecord(rAux);

					recBuf.head = headBuf.data();
					recBuf.headLen = rec.headLen;
					parser_2.WriteNextRecord(recBuf);
				}
				else
				{
					recBuf.head = rec.head;
					recBuf.headLen = rec.headLen;
					parser_1.WriteNextRecord(recBuf);

					rAux.head = headBuf.data();
					rAux.headLen = rec.headLen;
					parser_2.WriteNextRecord(rAux);
				}
			}
			else
			{
				if (rec.IsPairSwapped())
				{
					FastqRecord rAux = rec;
					rAux.head = headBuf.data();
					rAux.headLen = rec.headLen;

					parser_1.WriteNextRecord(rec.GetPair());
					parser_2.WriteNextRecord(rAux);
				}
				else
				{
					FastqRecord rAux = rec.GetPair();
					rAux.head = headBuf.data();
					rAux.headLen = rec.headLen;

					parser_1.WriteNextRecord(rec);
					parser_2.WriteNextRecord(rAux);
				}
			}

		}
		else
		{
			// use the same read name for both /1 and /2 reads
			//
			FastqRecord r = rec;
			r.head = headBuf.data();
			r.headLen = autoHeaderPrefix.size();
			r.headLen += to_string(headBuf.begin() + autoHeaderPrefix.size(), recStartIdx_);		// <library>.<num>

			if (r.IsReadReverse())
			{
				r.ComputeRC(recBuf);
				recBuf.head = r.head;
				recBuf.headLen = r.headLen;

				if (r.IsPairSwapped())
				{
					parser_1.WriteNextRecord(recBuf.GetPair());
					parser_2.WriteNextRecord(recBuf);
				}
				else
				{
					parser_1.WriteNextRecord(recBuf);
					parser_2.WriteNextRecord(recBuf.GetPair());
				}
			}
			else
			{
				if (r.IsPairSwapped())
				{
					parser_1.WriteNextRecord(r.GetPair());
					parser_2.WriteNextRecord(r);
				}
				else
				{
					parser_1.WriteNextRecord(r);
					parser_2.WriteNextRecord(r.GetPair());
				}
			}
		}

		recStartIdx_++;
	}

	chunk_.chunks[0]->size = parser_1.FinishParsing(SingleDnaRecordParser::ParseWrite);
	chunk_.chunks[1]->size = parser_2.FinishParsing(SingleDnaRecordParser::ParseWrite);
	ASSERT(chunk_.chunks[0]->size == chunk_.chunks[1]->size);

	return chunk_.chunks[0]->size;
}


uint64 FastqRecordsParserPE::ParseFrom(IFastqChunkCollection& chunk_, std::vector<FastqRecord>& records_, FastqRawBlockStats& stats_, bool keepComments_)
{
	ASSERT(chunk_.chunks.size() >= 3);

	FastqChunk* chunk1 = chunk_.chunks[FastqChunkCollectionPE::InputChunk1];
	FastqChunk* chunk2 = chunk_.chunks[FastqChunkCollectionPE::InputChunk2];
	FastqChunk* outChunk = chunk_.chunks[FastqChunkCollectionPE::OutputChunk];

	ASSERT(chunk1->size == chunk2->size);
	if (outChunk->data.Size() < chunk1->data.Size() + chunk2->data.Size())
		outChunk->data.Extend(chunk1->data.Size() + chunk2->data.Size());

	records_.clear();
#if EXTRA_MEM_OPT
	records_.shrink_to_fit();
#endif
	stats_.Clear();
	FastqRawBlockStats stats_2;

	SingleFastqRecordParser parser1(useHeaders, keepComments_), parser2(useHeaders, keepComments_);
	parser1.StartParsing(*chunk1, SingleDnaRecordParser::ParseRead, &stats_);
	parser2.StartParsing(*chunk2, SingleDnaRecordParser::ParseRead, &stats_2);

	FastqRecord rec1, rec2;
	char* outBufferPtr = (char*)outChunk->data.Pointer();

	while (parser1.ReadNextRecord(rec1) && parser2.ReadNextRecord(rec2))
	{
		ASSERT(rec1.seqLen == rec2.seqLen);

		// copy the sequence and update the pointer
		//
		std::copy(rec1.seq, rec1.seq + rec1.seqLen, outBufferPtr);
		rec1.seq = outBufferPtr;
		outBufferPtr += rec1.seqLen;

		std::copy(rec2.seq, rec2.seq + rec2.seqLen, outBufferPtr);
		outBufferPtr += rec2.seqLen;

		// copy the quality and update the pointer
		//
		std::copy(rec1.qua, rec1.qua + rec1.seqLen, outBufferPtr);
		rec1.qua = outBufferPtr;
		outBufferPtr += rec1.seqLen;

		std::copy(rec2.qua, rec2.qua + rec2.seqLen, outBufferPtr);
		outBufferPtr += rec2.seqLen;


		// store the record
		//
		rec1.auxLen = rec2.seqLen;
		records_.emplace_back(rec1);


		// simple validate headers (if used)
		//
		if (useHeaders)
		{
			ASSERT(stats_.head.fields.size() == stats_2.head.fields.size());
		}
	}
	ASSERT(records_.size() > 0);

	int64 bufSize = outBufferPtr - (char*)outChunk->data.Pointer();
	uint64 size_1 = parser1.FinishParsing(SingleDnaRecordParser::ParseRead);
	uint64 size_2 = parser2.FinishParsing(SingleDnaRecordParser::ParseRead);

	// this is an estimated size, as it does not take into account the skipped headers size
	//
	ASSERT(bufSize <= (int64)size_1 + (int64)size_2);
	outChunk->size = bufSize;


	// merge stats
	//
	stats_.minAuxLen = stats_2.minSeqLen;
	stats_.maxAuxLen = stats_2.maxSeqLen;

	// here we also update the stats regrding headers in PE mode, yet we'll select
	// the field coding the PE flag at the end of processing the archive
	stats_.Update(stats_2);


	return outChunk->size;
}


// dynamic memory management
//

// full FASTQ parsers
//
uint64 FastqRecordsParserDynSE::ParseTo(const std::vector<FastqRecord> &reads_, IFastqChunkCollection &chunk_, uint64 recStartIdx_)
{
	FastqRecordBuffer rcBuf;

	std::array<char, FastqRecord::MaxTagLen> headBuf;
	if (!useHeaders)
		std::copy(autoHeaderPrefix.begin(), autoHeaderPrefix.end(), headBuf.begin());


	// initialize chunks
	//
	if (chunk_.chunks.size() == 0)
	{
		chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));
	}
	else
	{
		for (DataChunk* c : chunk_.chunks)
		{
			if (c->data.Size() < DefaultBufferSize)
				c->data.Extend(DefaultBufferSize);
			c->size = 0;
		}
	}

	SingleFastqRecordParser parser(useHeaders);
	uint64 totalSize = 0;
	uint64 currentBufferPos = 0;
	uint64 currentChunkId = 0;

	DataChunk* currentChunk = chunk_.chunks[currentChunkId];
	parser.StartParsing(*currentChunk, SingleDnaRecordParser::ParseWrite);

	for (const auto& rec : reads_)
	{
		FastqRecord r;

		if (rec.IsReadReverse())
		{
			rec.ComputeRC(rcBuf);
			if (!useHeaders)
			{
				rcBuf.head = headBuf.data();
				rcBuf.headLen = autoHeaderPrefix.size();
				rcBuf.headLen += to_string(headBuf.begin() + autoHeaderPrefix.size(), recStartIdx_);		// <library>.<num>
			}
			else
			{
				rcBuf.head = rec.head;
				rcBuf.headLen = rec.headLen;
			}

			//parser.WriteNextRecord(rcBuf);
			r = rcBuf;
		}
		else
		{
			if (!useHeaders)
			{
				//FastqRecord r = rec;

				r = rec;
				r.head = headBuf.data();
				r.headLen = autoHeaderPrefix.size();
				r.headLen += to_string(headBuf.begin() + autoHeaderPrefix.size(), recStartIdx_);		// <library>.<num>

				//parser.WriteNextRecord(r);
			}
			else
			{
				//parser.WriteNextRecord(rec);
				r = rec;
			}
		}


		uint64 readSize = r.headLen + r.seqLen * 2 + 4 + 1;

		if (currentBufferPos + readSize > DefaultBufferSize)
		{
			currentChunk->size = parser.FinishParsing(SingleDnaRecordParser::ParseWrite);
			totalSize += currentChunk->size;

			currentChunkId++;

			if (chunk_.chunks.size() < currentChunkId + 1)
				chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));
			currentChunk = chunk_.chunks[currentChunkId];

			currentBufferPos = 0;

			parser.StartParsing(*currentChunk, SingleDnaRecordParser::ParseWrite);
		}

		parser.WriteNextRecord(r);
		currentBufferPos += readSize;

		recStartIdx_++;
	}

	//chunk_.chunks[0]->size = parser.FinishParsing(SingleDnaRecordParser::ParseWrite);

	if (currentBufferPos > 0)
	{
		currentChunk->size =  parser.FinishParsing(SingleDnaRecordParser::ParseWrite);
		totalSize += currentChunk->size;
	}


#if EXTRA_MEM_OPT

	std::vector<DataChunk*> validChunks;
	for (DataChunk* dc : chunk_.chunks)
	{
		if (dc->size > 0)
			validChunks.push_back(dc);
		else
			delete dc;
	}
	std::swap(chunk_.chunks, validChunks);

#endif

	return totalSize;
}


uint64 FastqRecordsParserDynPE::ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_)
{
	//ASSERT(chunk_.chunks.size() >= 2);

	SingleFastqRecordParser parser_1(useHeaders), parser_2(useHeaders);

	FastqRecordBuffer recBuf;

	std::array<char, FastqRecord::MaxTagLen> headBuf;

	if (!useHeaders)
		std::copy(autoHeaderPrefix.begin(), autoHeaderPrefix.end(), headBuf.begin());


	// initialize chunks
	//
	//ASSERT(chunk_.chunks.size() % 2 == 0);
	if (chunk_.chunks.size() % 2 == 1)
		chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));

	if (chunk_.chunks.size() == 0)
	{
		chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));
		chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));
	}
	else
	{
		for (DataChunk* c : chunk_.chunks)
		{
			if (c->data.Size() < DefaultBufferSize)
				c->data.Extend(DefaultBufferSize);
			c->size = 0;
		}
	}



	uint64 totalSize_1 = 0, totalSize_2 = 0;
	uint64 currentBufferPos_1 = 0, currentBufferPos_2 = 0;
	uint64 currentChunkId_1 = 0, currentChunkId_2 = 1;

	DataChunk* currentChunk_1 = chunk_.chunks[currentChunkId_1];
	DataChunk* currentChunk_2 = chunk_.chunks[currentChunkId_2];

	parser_1.StartParsing(*currentChunk_1, SingleDnaRecordParser::ParseWrite);
	parser_2.StartParsing(*currentChunk_2, SingleDnaRecordParser::ParseWrite);


	for (const FastqRecord& rec : reads_)
	{
		FastqRecord r1, r2;


		if (useHeaders)
		{
			// here we set the proper headers
			//
			std::copy(rec.head, rec.head + rec.headLen, headBuf.begin());


			// check whether the /1 /2 fields are present
			//
			if (peFieldIdx != 0)
			{
				// find position of the /2 token
				const std::string separators = FastqRawBlockStats::HeaderStats::Separators();
				ASSERT(peFieldIdx != 0);
				uint32 pairTokenPos = 0;

				uint32 fidx = 0;
				for (uint32 i = 0; i <= rec.headLen; ++i)
				{
					if (!std::count(separators.begin(), separators.end(), rec.head[i]) && (i != rec.headLen))
						continue;

					fidx++;

					if (fidx == peFieldIdx)
					{
						pairTokenPos = i+1;
						break;
					}
				}

				ASSERT(rec.head[pairTokenPos] == '1');
				headBuf[pairTokenPos] = '2';
			}


			if (rec.IsReadReverse())
			{
				rec.ComputeRC(recBuf);
				FastqRecord rAux = recBuf.GetPair();

				if (rec.IsPairSwapped())
				{
					rAux.head = rec.head;
					rAux.headLen = rec.headLen;
					//parser_1.WriteNextRecord(rAux);
					r1 = rAux;

					recBuf.head = headBuf.data();
					recBuf.headLen = rec.headLen;
					//parser_2.WriteNextRecord(recBuf);
					r2 = recBuf;
				}
				else
				{
					recBuf.head = rec.head;
					recBuf.headLen = rec.headLen;
					//parser_1.WriteNextRecord(recBuf);
					r1 = recBuf;

					rAux.head = headBuf.data();
					rAux.headLen = rec.headLen;
					//parser_2.WriteNextRecord(rAux);
					r2 = rAux;
				}
			}
			else
			{
				if (rec.IsPairSwapped())
				{
					FastqRecord rAux = rec;
					rAux.head = headBuf.data();
					rAux.headLen = rec.headLen;

					//parser_1.WriteNextRecord(rec.GetPair());
					//parser_2.WriteNextRecord(rAux);

					r1 = rec.GetPair();
					r2 = rAux;
				}
				else
				{
					FastqRecord rAux = rec.GetPair();
					rAux.head = headBuf.data();
					rAux.headLen = rec.headLen;

					//parser_1.WriteNextRecord(rec);
					//parser_2.WriteNextRecord(rAux);
					r1 = rec;
					r2 = rAux;
				}
			}

		}
		else
		{
			// use the same read name for both /1 and /2 reads
			//
			FastqRecord r = rec;
			r.head = headBuf.data();
			r.headLen = autoHeaderPrefix.size();
			r.headLen += to_string(headBuf.begin() + autoHeaderPrefix.size(), recStartIdx_);		// <library>.<num>

			if (r.IsReadReverse())
			{
				r.ComputeRC(recBuf);
				recBuf.head = r.head;
				recBuf.headLen = r.headLen;

				if (r.IsPairSwapped())
				{
					//parser_1.WriteNextRecord(recBuf.GetPair());
					//parser_2.WriteNextRecord(recBuf);
					r1 = recBuf.GetPair();
					r2 = recBuf;
				}
				else
				{
					//parser_1.WriteNextRecord(recBuf);
					//parser_2.WriteNextRecord(recBuf.GetPair());
					r1 = recBuf;
					r2 = recBuf.GetPair();
				}
			}
			else
			{
				if (r.IsPairSwapped())
				{
					//parser_1.WriteNextRecord(r.GetPair());
					//parser_2.WriteNextRecord(r);
					r1 = r.GetPair();
					r2 = r;
				}
				else
				{
					//parser_1.WriteNextRecord(r);
					//parser_2.WriteNextRecord(r.GetPair());
					r1 = r;
					r2 = r.GetPair();
				}
			}
		}



		// store r_1
		uint64 readSize_1 = r1.headLen + r1.seqLen * 2 + 4 + 1;

		if (currentBufferPos_1 + readSize_1 > DefaultBufferSize)
		{
			currentChunk_1->size = parser_1.FinishParsing(SingleDnaRecordParser::ParseWrite);
			totalSize_1 += currentChunk_1->size;

			// we need to increment the IDs in pairs, not to overlap the buffers
			currentChunkId_1 += 2;

			while (chunk_.chunks.size() < currentChunkId_1 + 1 + 1)	// extra +1 for PE buffer
				chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));
			ASSERT(chunk_.chunks.size() % 2 == 0);

			currentChunk_1 = chunk_.chunks[currentChunkId_1];

			currentBufferPos_1 = 0;

			parser_1.StartParsing(*currentChunk_1, SingleDnaRecordParser::ParseWrite);
		}

		parser_1.WriteNextRecord(r1);
		currentBufferPos_1 += readSize_1;


		// store r_2
		//

		uint64 readSize_2 = r2.headLen + r2.seqLen * 2 + 4 + 1;
		if (currentBufferPos_2 + readSize_2 > DefaultBufferSize)
		{
			currentChunk_2->size = parser_2.FinishParsing(SingleDnaRecordParser::ParseWrite);
			totalSize_2 += currentChunk_2->size;

			// we need to increment Ids in pairs
			currentChunkId_2 += 2;

			while (chunk_.chunks.size() < currentChunkId_2 + 1)
				chunk_.chunks.push_back(new DataChunk(DefaultBufferSize));
			ASSERT(chunk_.chunks.size() % 2 == 0);

			currentChunk_2 = chunk_.chunks[currentChunkId_2];

			currentBufferPos_2 = 0;

			parser_2.StartParsing(*currentChunk_2, SingleDnaRecordParser::ParseWrite);
		}

		parser_2.WriteNextRecord(r2);
		currentBufferPos_2 += readSize_2;


		recStartIdx_++;
	}

	if (currentBufferPos_1 > 0)
	{
		currentChunk_1->size = parser_1.FinishParsing(SingleDnaRecordParser::ParseWrite);
		totalSize_1 += currentChunk_1->size;
	}

	if (currentBufferPos_2 > 0)
	{
		currentChunk_2->size = parser_2.FinishParsing(SingleDnaRecordParser::ParseWrite);
		totalSize_2 += currentChunk_2->size;
	}

	ASSERT(totalSize_1 == totalSize_2);



#if EXTRA_MEM_OPT

	std::vector<DataChunk*> validChunks;
	for (size_t i = 0; i < chunk_.chunks.size(); i += 2)
	{
		ASSERT(chunk_.chunks.size() > i + 1);

		DataChunk* c1 = chunk_.chunks[i];
		DataChunk* c2 = chunk_.chunks[i+1];

		if (c1->size == 0 && c2->size == 0)
		{
			delete c1;
			delete c2;
		}
		else
		{
			validChunks.push_back(c1);
			validChunks.push_back(c2);
		}
	}
	std::swap(chunk_.chunks, validChunks);

#endif

	return totalSize_1;
}
