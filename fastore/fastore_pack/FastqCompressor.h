/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM and QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINCOMPRESSOR
#define H_BINCOMPRESSOR

#include "../fastore_bin/Globals.h"

#include <string>
#include <deque>
#include <vector>
#include <queue>
#include <unordered_map>

#include "Params.h"
#include "CompressedBlockData.h"
#include "ReadsClassifier.h"
#include "ContigBuilder.h"

#include "../fastore_bin/BitMemory.h"
#include "../fastore_bin/Params.h"
#include "../fastore_bin/FastqRecord.h"
#include "../fastore_bin/FastqCategorizer.h"
#include "../rle/RleEncoder.h"
#include "../rc/ContextEncoder.h"
#include "../ppmd/PPMd.h"

#include "../fastore_rebin/Params.h"
#include "qv_compressor.h"

#include "../fastore_bin/QVZ.h"


#define ENC_HR_AC 0


/**
 * Basic compressor/decompressor interfaces
 */
class IStoreBase
{
public:
	IStoreBase(const CompressorParams& params_, const CompressorAuxParams& auxParams_);
	virtual ~IStoreBase();

protected:
	struct BaseBlockHeader
	{
		static const uint64 Size = 4*sizeof(uint64) + 2*sizeof(uint32) + 2*sizeof(uint8);

		uint64 recordsCount;
		uint64 rawDnaStreamSize;
		uint64 rawIdStreamSize;
		uint64 footerOffset;

		uint32 footerSize;
		uint32 minimizerId;

		uint8 recMinLen;
		uint8 recMaxLen;


		BaseBlockHeader()
		{
			Reset();
		}

		void Reset()
		{
			minimizerId = 0;
			recordsCount = 0;

			recMinLen = (uint8)-1;
			recMaxLen = 0;

			rawDnaStreamSize = 0;
			rawIdStreamSize = 0;

			footerSize = footerOffset = 0;
		}
	};

	struct BaseBlockFooter
	{
		// extra quality scores compression information
		//
		uint32 sampleValue;		// some dummy value placeholder


		BaseBlockFooter()
			:	sampleValue(0)
		{}

		void Reset()
		{
			sampleValue = 0;
		}
	};


	static const uint64 DefaultPpmdMemorySizeMb = 16;
	static const uint32 DefaultPpmdOrder = 4;

	const CompressorParams params;				// TODO: try ref
	const CompressorAuxParams auxParams;


	// for dry run
	Buffer* dryFastqBuffer;
	BitMemoryWriter* dryFastqWriter;
	uint64 currentRecordIdx;


	void StoreRawHeader(const BaseBlockHeader& header_, BitMemoryWriter& writer_);
	void ReadRawHeader(BaseBlockHeader& header_, BitMemoryReader& reader_);

	void StoreRawFooter(const BaseBlockFooter& footer_, BitMemoryWriter& writer_);
	void ReadRawFooter(BaseBlockFooter& footer_, BitMemoryReader& reader_);

	void CompressBuffer(PpmdEncoder& encoder_, const DataChunk& inChunk_, uint64 inSize_,
						DataChunk& outChunk_, uint64& outSize_, uint64 outOffset_);

	void DecompressBuffer(PpmdDecoder& decoder_, DataChunk& outChunk_, uint64& outSize_,
						  const DataChunk& inChunk_, uint64 inSize_, uint64 inOffset_);
};


class IDnaStoreBase
{
public:
	IDnaStoreBase(const MinimizerParameters &minimizer_);

protected:
	std::array<char, 128> dnaToIdx;
	std::array<char, 8> idxToDna;
};


class IQualityStoreBase
{
public:
	IQualityStoreBase(const QualityCompressionData& globalQuaData_);

protected:
	// stats for quality compression
	//
	const QualityCompressionData& globalQuaData;
	well_state_t local_well;


	// default quality coders
	//
	typedef TAdvancedContextCoder<2, 10> BinaryQuaCoder;
	typedef TAdvancedContextCoder<8, 6> Illu8QuaCoder;

	struct QualityEncoders
	{
		typedef TEncoder<BinaryQuaCoder> BinaryEncoder;
		typedef TEncoder<Illu8QuaCoder> Illu8Encoder;
		//typedef TEncoder<QVZQuaCoder> QVZEncoder;

		BinaryEncoder* binaryCoder;
		Illu8Encoder* illu8Coder;
		BitMemoryWriter* rawCoder;
		QVZEncoder* qvzCoder;

		QualityEncoders()
			:	binaryCoder(NULL)
			,	illu8Coder(NULL)
			,	rawCoder(NULL)
			,   qvzCoder(NULL)
		{}
	};

	struct QualityDecoders
	{
		typedef TDecoder<BinaryQuaCoder> BinaryDecoder;
		typedef TDecoder<Illu8QuaCoder> Illu8Decoder;
		//typedef TDecoder<QVZQuaCoder> QVZEncoder;

		BinaryDecoder* binaryCoder;
		Illu8Decoder* illu8Coder;
		BitMemoryReader* rawCoder;
		QVZDecoder* qvzCoder;

		QualityDecoders()
			:	binaryCoder(NULL)
			,	illu8Coder(NULL)
			,	rawCoder(NULL)
			,   qvzCoder(NULL)
		{}
	};

	std::array<char, 64> quaToIdx_8bin;
	std::array<char, 8> idxToQua_8bin;


	// INFO: a robust way to handle quality compression in one place
	// assuming that the same/similar technique is (so far) used
	// to compress reads in standard, small and N' bins
	void CompressReadQuality(const FastqRecord& rec_,
							 QualityEncoders& enc_,
							 const CompressorParams& params_,
							 const bool dryRun_ = false,
							 std::vector<byte>* dryBuffer_ = NULL);

	void DecompressReadQuality(FastqRecord& rec_,
							   QualityDecoders& dec_,
							   const CompressorParams& params_);

	void ResetWellRng();
};


class IHeaderStoreBase
{
public:
	IHeaderStoreBase(const FastqRawBlockStats::HeaderStats& headData_)
		:	headData(headData_)
	{
		SetupFieldCompressionSpec();
	}

protected:
	struct FieldCompressionSpec
	{
		enum CompressionMethod
		{
			COMP_CONST,		// const value, no need for compression
			COMP_TOKEN,		// selection of token from range - can use Huffman or range-coder
			COMP_RAW		// store raw numeric value
		};

		CompressionMethod method;
		std::vector<std::string> tokenValues;	// a copy for faster access
		uint8 bitsPerValue;

		FieldCompressionSpec()
			:	method(COMP_RAW)
			,	bitsPerValue(0)
		{}
	};

	typedef TAdvancedContextCoder<256, 1> TokenCoder;
	typedef TAdvancedContextCoder<256, 1> ValueCoder;

	struct FieldEncoders
	{
		typedef TEncoder<TokenCoder> TokenEncoder;
		typedef TEncoder<ValueCoder> ValueEncoder;

		TokenEncoder* tokenCoder;
		ValueEncoder* valueCoder;

		FieldEncoders()
			:	tokenCoder(NULL)
			,	valueCoder(NULL)
		{}
	};

	struct FieldDecoders
	{
		typedef TDecoder<TokenCoder> TokenDecoder;
		typedef TDecoder<ValueCoder> ValueDecoder;

		TokenDecoder* tokenCoder;
		ValueDecoder* valueCoder;

		FieldDecoders()
			:	tokenCoder(NULL)
			,	valueCoder(NULL)
		{}
	};

	const FastqRawBlockStats::HeaderStats& headData;

	std::vector<FieldCompressionSpec> fieldsSpec;

	void SetupFieldCompressionSpec();

	void CompressReadId(const FastqRecord& rec_,
						FieldEncoders& enc_) const;

	void DecompressReadId(FastqRecord& rec_,
						  FieldDecoders& dec_) const;
};


/**
 * Basic functionality to store the reads in raw form,
 * i.e. not differentially compressed ('N' bin)
 *
 */
class IDnaRawStoreBase :	public IStoreBase,
							public IDnaStoreBase,
							public IQualityStoreBase,
							public IHeaderStoreBase
{
public:
	IDnaRawStoreBase(const CompressorParams& params_,
					 const QualityCompressionData& globalQuaData_,
					 const FastqRawBlockStats::HeaderStats& headData_,
					 const CompressorAuxParams& auxParams_ = CompressorAuxParams())
		:	IStoreBase(params_, auxParams_)
		,	IDnaStoreBase(params_.minimizer)
		,	IQualityStoreBase(globalQuaData_)
		,	IHeaderStoreBase(headData_)
	{}

protected:
	struct RawBlockHeader : public BaseBlockHeader
	{
		static const uint64 Size = BaseBlockHeader::Size + 4*sizeof(uint64);

		// TODO: can be added more buffers as in LzBlockHeader
		// and used as only one unified structure
		uint64 dnaCompSize;
		uint64 quaCompSize;
		uint64 idTokenCompSize;
		uint64 idValueCompSize;

		RawBlockHeader()
			:	dnaCompSize(0)
			,	quaCompSize(0)
			,	idTokenCompSize(0)
			,	idValueCompSize(0)
		{}

		void Reset()
		{
			BaseBlockHeader::Reset();

			dnaCompSize = 0;
			quaCompSize = 0;
			idTokenCompSize = 0;
			idValueCompSize = 0;
		}
	};


	struct CompressedBlockDescription
	{
		RawBlockHeader header;
		BaseBlockFooter footer;
		bool isLenConst;

		CompressedBlockDescription()
		{
			Reset();
		}

		void Reset()
		{
			header.Reset();
			footer.Reset();
			isLenConst = true;
		}
	};

	CompressedBlockDescription blockDesc;


	void StoreHeader(const RawBlockHeader& header_, DataChunk& buffer_);
	void ReadHeader(RawBlockHeader& header_, DataChunk& buffer_);
};


/**
 * Basic functionality to store the reads encoded differentially
 *
 */
class ILzCompressorBase :	public IStoreBase,
							public IDnaStoreBase,
							public IQualityStoreBase,
							public IHeaderStoreBase
{
public:
	ILzCompressorBase(const CompressorParams& params_,
					  const QualityCompressionData& globalQuaData_,
					  const FastqRawBlockStats::HeaderStats& headData_,
					  const CompressorAuxParams& auxParams_ = CompressorAuxParams())
		:	IStoreBase(params_, auxParams_)
		,	IDnaStoreBase(params_.minimizer)
		,	IQualityStoreBase(globalQuaData_)
		,	IHeaderStoreBase(headData_)
	{}

protected:
	enum ReadFlags
	{
		ReadIdentical = 0,
		ReadDifficult,

		ReadShiftOnly,
		ReadFullEncode,
		ReadFullExpensive,

		ReadTreeGroupStart,

		ReadContigGroupStart,
		ReadContigGroupNext,

		ReadGroupEnd
	};

	enum ReadMatchType
	{
		HardRead,
		ExactMatch,
		LzMatch,
		ContigRead
	};


	struct LzBlockHeader : public BaseBlockHeader
	{
		std::vector<uint64> workBufferSizes;
		std::vector<uint64> compBufferSizes;

		void Reset(uint32 buffersCount_)
		{
			BaseBlockHeader::Reset();

			workBufferSizes.resize(buffersCount_);
			compBufferSizes.resize(buffersCount_);

			std::fill(workBufferSizes.begin(), workBufferSizes.end(), 0);
			std::fill(compBufferSizes.begin(), compBufferSizes.end(), 0);
		}

		static uint64 Size(uint64 buffersNum_)
		{
			return BaseBlockHeader::Size + 2 * buffersNum_ * sizeof(uint64);
		}
	};

	struct CompressedBlockDescription
	{
		LzBlockHeader header;
		BaseBlockFooter footer;
		bool isLenConst;

		CompressedBlockDescription()
			:	isLenConst(true)
		{}

		void Reset(uint32 buffersCount_)
		{
			header.Reset(buffersCount_);
			footer.Reset();
			isLenConst = true;
		}
	};

	struct LzContext
	{
		std::array<char, MAX_SIGNATURE_LEN> currentMinimizerBuf;

		std::vector<const FastqRecord*> history;
		uint32 lastBufPos;

		LzContext()
			:	lastBufPos(0)
		{
			std::fill(currentMinimizerBuf.begin(), currentMinimizerBuf.end(), 0);
		}
	};


	// dna coders
	//
	typedef TSimpleContextCoder<2, 4> RevContextCoder;
	typedef TAdvancedContextCoder<8, 4> LettersContextCoder;
	typedef TSimpleContextCoder<8, 4> HardContextCoder;

	struct DnaEncoders
	{
		typedef TEncoder<RevContextCoder> RevEncoder;
		typedef TEncoder<LettersContextCoder> LettersEncoder;
		typedef TEncoder<HardContextCoder> HardEncoder;

		BinaryRleEncoder* matchRleCoder;
		BinaryRleEncoder* consMatchCoder;
		Rle0Encoder* lzRle0Coder;

		//FlagEncoder* flagCoder;
		RevEncoder* revRcCoder;
		RevEncoder* matchRcCoder;

		LettersEncoder* lettersXRcCoder;
		LettersEncoder* lettersCRcCoder;
		HardEncoder* hardCoder;

		DnaEncoders()
			:	matchRleCoder(NULL)
			,	lzRle0Coder(NULL)
		//	,	flagCoder(NULL)
			,	revRcCoder(NULL)
			,	matchRcCoder(NULL)
			,	lettersXRcCoder(NULL)
			,	lettersCRcCoder(NULL)
			,	hardCoder(NULL)
		{}
	};

	struct DnaDecoders
	{
		//typedef TDecoder<FlagContextCoder> FlagDecoder;
		typedef TDecoder<RevContextCoder> RevDecoder;
		typedef TDecoder<LettersContextCoder> LettersDecoder;
		typedef TDecoder<HardContextCoder> HardDecoder;

		BinaryRleDecoder* matchRleCoder;
		BinaryRleDecoder* consMatchCoder;
		Rle0Decoder* lzRle0Coder;

		//FlagDecoder* flagCoder;
		RevDecoder* revRcCoder;
		RevDecoder* matchRcCoder;
		LettersDecoder* lettersXRcCoder;
		LettersDecoder* lettersCRcCoder;
		HardDecoder* hardCoder;

		DnaDecoders()
			:	matchRleCoder(NULL)
			,	consMatchCoder(NULL)
			,	lzRle0Coder(NULL)
			//,	flagCoder(NULL)
			,	revRcCoder(NULL)
			,	matchRcCoder(NULL)
			,	lettersXRcCoder(NULL)
			,	lettersCRcCoder(NULL)
			,	hardCoder(NULL)
		{}
	};


	// predefined constants
	//
	static const int32 ShiftOffset = 128 + 1;
	static const uint32 InitialSeqLen = 255;				// TODO: max_seq_len -- define
	static const uint32 MinimizerPositionSymbol = '.';


	CompressedBlockDescription blockDesc;
	std::stack<LzContext> lzContextStack;

	void StoreHeader(const LzBlockHeader& header_, DataChunk& buffer_);
	void ReadHeader(LzBlockHeader& header_, DataChunk& buffer_);

	virtual void SetupBufferMask(std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_)
	{
		// setup decoding mask -- TODO: create only once and reuse
		//
		ppmdBufferCompMask_.resize(buffersNum_);
		std::fill(ppmdBufferCompMask_.begin(), ppmdBufferCompMask_.end(), true);
		//bufferMask[FastqWorkBuffersSE::FlagBuffer] = false;
		ppmdBufferCompMask_[FastqWorkBuffersSE::RevBuffer] = false;
		ppmdBufferCompMask_[FastqWorkBuffersSE::MatchBinaryBuffer] = false;
		ppmdBufferCompMask_[FastqWorkBuffersSE::LetterXBuffer] = false;
		ppmdBufferCompMask_[FastqWorkBuffersSE::ConsensusLetterBuffer] = false;
#if (ENC_HR_AC)
		ppmdBufferCompMask_[FastqWorkBuffersSE::HardReadsBuffer] = false;
#endif

		// WARN: in any other case than lossless compression, quality will
		// be compressed using own coder
		//
		if (params.quality.method != QualityCompressionParams::MET_NONE)
			ppmdBufferCompMask_[FastqWorkBuffersSE::QualityBuffer] = false;

		// do not compress the read id buffer
		ppmdBufferCompMask_[FastqWorkBuffersSE::ReadIdTokenBuffer] = false;
		ppmdBufferCompMask_[FastqWorkBuffersSE::ReadIdValueBuffer] = false;
	}
};



/**
 * Single-end FASTQ reads compressors/decompressors
 *
 */
class LzCompressorSE : public ILzCompressorBase
{
public:
	LzCompressorSE(const CompressorParams& params_,
				   const QualityCompressionData& globalQuaData_,
				   const FastqRawBlockStats::HeaderStats& headData_,
				   const CompressorAuxParams& auxParams_
				   /*, PPMdCompressor& ...*/);

	~LzCompressorSE();

	virtual void Compress(const std::vector<FastqRecord>& records_,
							 PackContext& packCtx_,
							 uint32 minimizerId_,
							 uint64 rawDnaStreamSize_,
							 FastqCompressedBin& dnaWorkBin_,
							 CompressedFastqBlock &compBin_);

protected:
	struct ConsensusEncoder
	{
		enum GroupType
		{
			ConsensusGroup = 0,
			TreeGroup
		};

		GroupType type;
		int32 lastMinimPos;
		const ConsensusDefinition* definition;

		ConsensusEncoder()
			:	type(ConsensusGroup)
			,	lastMinimPos(0)
			,	definition(NULL)
		{}
	};


	struct EncodingContext
	{
		DnaEncoders dna;
		QualityEncoders qua;
		FieldEncoders id;

		std::vector<BitMemoryWriter*> writers;
		std::vector<bool> bufferCompMask;
	};


	//typedef TEncoder<FlagContextCoder> FlagEncoder;

	static const int32 MaxShiftValue = ShiftOffset-1;


	// methods to be overloaded by PE
	//
	virtual void CompressMatch(const MatchNode* node_, LzContext& lzEncoder_);
	virtual void CompressRead(const FastqRecord& record_, ReadMatchType type_, bool aux_ = false);

	virtual void DryRunCompress(const FastqRecord& record_, BitMemoryWriter* dryFastqWriter_);


	virtual void StartEncoding(std::vector<DataChunk*>& buffers_);
	virtual void EndEncoding(std::vector<DataChunk*>& buffers_);

	virtual void SetupBuffers(Buffer& buffer_, std::vector<bool>& ppmdBufferCompMask_, uint32 bufferNum_);

	void CompressBuffers(std::vector<DataChunk*>& buffers_,
						 const std::vector<bool>& ppmdBufferCompMask_,		// TODO: make constant in ctor
						 DataChunk &compBin_,
						 uint64 chunkOffset_);

	void CompressRecords(PackContext& packCtx_);


	// experimental quality compression methods
	//
	// TODO: ideally can be implemented as an external set of module
	//
	virtual void CompressReadQuality(const FastqRecord& record_);

	virtual void CompressQuality();



	EncodingContext mainCtx;

	ReadsClassifierSE readsClassifier;


	// tmp
	CompressedFastqBlockStats* blockStats;
	//

private:
	void CompressNode(MatchNode* node_);


	// single records compression
	//
	void CompressNormalMatch(const MatchNode* node_, uint32 lzId_, bool expensiveEncoding_ = false);

	void CompressExactRead(const FastqRecord& record_);

	void CompressHardRead(const FastqRecord& record_);


	// consensus compression
	//
	void CompressContig(MatchNode* node_);

	void StoreContigDefinition(const ConsensusDefinition& contig_,
							   const char* mainSeq_,
							   uint32 mainSignaturePos_);

	void CompressContigRecords(MatchNode* mainNode_, const ContigDefinition& contig_);

	void CompressContigRead(const FastqRecord& rec_, bool useTreeShiftBuffer_ = false);


	// other groups compression
	//
	void CompressSubTree(const MatchNode* rootNode_);

	void CompressExactChildren(const MatchNode* node_);


	PpmdEncoder* ppmdCoder;

	ContigBuilder consBuilder;

	std::stack<ConsensusEncoder> consEncoders;			// TODO: tree encoders, cons encoders
};


class LzDecompressorSE : public ILzCompressorBase
{
public:
	LzDecompressorSE(const CompressorParams& params_,
					 const QualityCompressionData& globalQuaData_,
					 const FastqRawBlockStats::HeaderStats& headData_);
	~LzDecompressorSE();

	virtual void Decompress(CompressedFastqBlock& compBin_,
							   std::vector<FastqRecord>& reads_,
							   FastqCompressedBin& fastqWorkBin_,
							   FastqChunk& dnaBuffer_);

protected:
	struct DecodingContext
	{
		DnaDecoders dna;
		QualityDecoders qua;
		FieldDecoders id;

		std::vector<BitMemoryReader*> readers;
		std::vector<bool> bufferCompMask;
	};


	virtual void DecompressRecord(FastqRecord& rec_, ReadFlags flag_,
								  const FastqRecord* auxRec_ = NULL, bool aux_ = false);
	virtual void DecompressMetaAndLinkRecords(std::vector<FastqRecord>& reads_, FastqChunk& chunk_);

	virtual void SetupBuffers(Buffer& outBuffer_, std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_);

	virtual void StartDecoding(std::vector<DataChunk*>& buffers_);
	virtual void EndDecoding();

	void DecompressBuffers(std::vector<DataChunk*>& buffers_,
						   const std::vector<bool>& ppmdBufferCompMask_,		// TODO: make constant or bind with buffers
						   DataChunk& compChunk_,
						   uint64 chunkOffset_);


	void DecompressRecords(std::vector<FastqRecord>& reads_);

	// experimental quality compression methods
	//
	// TODO: ideally can be implemented as an external set of module
	//
	virtual void DecompressReadQuality(FastqRecord& record_);
	virtual void DecompressQuality();


	DecodingContext mainCtx;

private:
	struct ConsensusDecoder
	{
		enum GroupType
		{
			ConsensusGroup,
			TreeGroup
		};

		GroupType type;
		FastqRecord* mainLzRecord;

		uint32 seqLen;					// deprecated

		uint32 readsCount;
		int32 lastMinimPos;				// used for delta decoding

		ConsensusDefinition definition;
		std::array<char, MAX_SIGNATURE_LEN> curSignatureString;

		ConsensusDecoder(uint32 seqLen_)
		{
			Reset(seqLen_);
		}

		void Reset(uint32 seqLen_)
		{
			type = ConsensusGroup;
			mainLzRecord = NULL;

			lastMinimPos = 0;
			readsCount = 0;

			if (seqLen != seqLen_ || definition.sequence.size() != seqLen_ * 2)
			{
				definition.sequence.resize(seqLen_ * 2);
				definition.variantPositions.resize(seqLen_ * 2);

				seqLen = seqLen_;
			}
			std::fill(definition.sequence.begin(), definition.sequence.end(), 'N');
			std::fill(definition.variantPositions.begin(), definition.variantPositions.end(), false);

			definition.range.first = seqLen;
			definition.range.second = seqLen;
		}
	};

	void DecompressExactRead(FastqRecord& rec_, const FastqRecord& aux_);
	void DecompressHardRead(FastqRecord& rec_);
	void DecompressLzMatch(FastqRecord& rec_, uint32 flag_);

	void DecompressConsensus(std::vector<FastqRecord>& reads_, uint64& startIdx_);		// TODO: here we should use iterator
	void DecompressConsensusRead(FastqRecord& rec_, bool useTreeShiftBuffer_);
	void DecompressConsensusExactMatch(FastqRecord& rec_, const FastqRecord& prev_);

	void DecompressTree(const FastqRecord& rootRec_,									// TODO: here we should use iterator
						std::vector<FastqRecord>& reads_,
						uint64& startIdx_);

	PpmdDecoder *ppmdCoder;

	std::stack<ConsensusDecoder> consDecoders;


	// DBG:
	uint32 currentSignature;
};

#if 0
class LzDecompressorDynSE : public LzDecompressorSE
{
public:
	using LzDecompressorSE::LzDecompressorSE;
	// TODO: initialize members


	virtual void StartDecompress(CompressedFastqBlock& compBin_,
								 FastqCompressedBin& fastqWorkBin_);

	virtual uint64 DecompressNext(std::vector<FastqRecord>& reads_,
								  FastqChunk& dnaBuffer_);

	virtual void FinishDecompress();

private:
	uint64 decompressedRecs;
	uint64 decompressedBytes;

	void SetupBuffers(Buffer& outBuffer_, std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_);
};
#endif

class RawCompressorSE : public IDnaRawStoreBase
{
public:
	RawCompressorSE(const CompressorParams& params_,
					const QualityCompressionData& globalQuaData_,
					const FastqRawBlockStats::HeaderStats& headData_,
					const CompressorAuxParams& auxParams_);
	~RawCompressorSE();

	virtual void Compress(const std::vector<FastqRecord>& records_,
							 PackContext& packCtx_,
							 uint32 minimizerId_,
							 uint64 rawDnaStreamSize_,
							 FastqCompressedBin& dnaWorkBin_,
							 CompressedFastqBlock &compBin_);

protected:
	struct EncodingContext
	{
		struct DnaEncoders
		{
			BitMemoryWriter* rawCoder;

			DnaEncoders()
				:	rawCoder(NULL)
			{}
		};

		BitMemoryWriter* dnaWriter;
		BitMemoryWriter* quaWriter;
		BitMemoryWriter* idTokenWriter;
		BitMemoryWriter* idValueWriter;

		DnaEncoders dna;
		QualityEncoders qua;
		FieldEncoders id;

		EncodingContext()
			:	dnaWriter(NULL)
			,	quaWriter(NULL)
			,	idTokenWriter(NULL)
			,	idValueWriter(NULL)
		{}
	};

	virtual void StartEncoding(std::vector<DataChunk*>& buffers_);
	virtual void EndEncoding(std::vector<DataChunk*>& buffers_);

	virtual void CompressReadSequence(const FastqRecord& rec_);
	virtual void CompressReadQuality(const FastqRecord& rec_);

	virtual void CompressDna(const DataChunk& inChunk_, uint64 inSize_,
							 DataChunk& outChunk_, uint64& outSize_, uint64 outOffset_);

	virtual void CompressQuality(const DataChunk& inChunk_, uint64 inSize_,
								 DataChunk& outChunk_, uint64& outSize_, uint64 outOffset_);

	virtual void DryRunCompress(const FastqRecord& record_, BitMemoryWriter* dryFastqWriter_);


	EncodingContext mainCtx;

	PpmdEncoder* ppmdCoder;
};


class RawDecompressorSE : public IDnaRawStoreBase
{
public:
	RawDecompressorSE(const CompressorParams& params_,
					  const QualityCompressionData& globalQuaData_,
					  const FastqRawBlockStats::HeaderStats& headData_);
	~RawDecompressorSE();

	virtual void Decompress(CompressedFastqBlock& compBin_,
							   std::vector<FastqRecord>& reads_,
							   FastqCompressedBin& fastqWorkBin_,
							   FastqChunk& dnaBuffer_);

protected:
	struct DecodingContext
	{
		struct DnaEncoders
		{
			BitMemoryReader* rawCoder;

			DnaEncoders()
				:	rawCoder(NULL)
			{}
		};

		BitMemoryReader* dnaReader;
		BitMemoryReader* quaReader;
		BitMemoryReader* idTokenReader;
		BitMemoryReader* idValueReader;

		DnaEncoders dna;
		QualityDecoders qua;
		FieldDecoders id;

		DecodingContext()
			:	dnaReader(NULL)
			,	quaReader(NULL)
			,	idTokenReader(NULL)
			,	idValueReader(NULL)
		{}
	};

	void StartDecoding(std::vector<DataChunk*>& buffers_);
	void EndDecoding();

	virtual void DecompressRecords(std::vector<FastqRecord>& reads_, FastqChunk& Buffer_);

	virtual void DecompressReadSequence(FastqRecord& rec_);
	virtual void DecompressReadQuality(FastqRecord& rec_);

	virtual void DecompressDna(DataChunk& outChunk_, uint64& outSize_,
							   const DataChunk& inChunk_, uint64 inSize_, uint64 inOffset_);

	virtual void DecompressQuality(DataChunk& outChunk_, uint64& outSize_,
								   const DataChunk& inChunk_, uint64 inSize_, uint64 inOffset_);


	DecodingContext mainCtx;

	PpmdDecoder* ppmdCoder;
};


/**
 * Paired-end FASTQ reads compressors/decompressors
 *
 */
class LzCompressorPE : public LzCompressorSE
{
public:
	LzCompressorPE(const CompressorParams& params_,
				   const QualityCompressionData& globalQuaData_,
				   const FastqRawBlockStats::HeaderStats& headData_,
				   const CompressorAuxParams& auxParams_);
	~LzCompressorPE();

	void Compress(const std::vector<FastqRecord>& records_,
								 PackContext& packCtx_,
								 uint32 minimizerId_,
								 uint64 rawDnaStreamSize_,
								 FastqCompressedBin& dnaWorkBin_,
								 CompressedFastqBlock &compBin_)
	{
		LzCompressorSE::Compress(records_, packCtx_, minimizerId_,
								 rawDnaStreamSize_, dnaWorkBin_, compBin_);

		ClearPairBuffer();
	}

private:
	enum ReadFlagsPE
	{
		ReadDifficultPE,
		ReadIdenticalPE,
		ReadFullEncodePE,
		ReadFullExpensivePE,
		ReadShiftOnlyPE,
	};

	void CompressMatch(const MatchNode* node_, LzContext& lzEncoder_);
	void CompressRead(const FastqRecord& record_, ReadMatchType type_, bool aux_ = false);

	void SetupBufferMask(std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_);
	void StartEncoding(std::vector<DataChunk*>& buffers_);
	void EndEncoding(std::vector<DataChunk*>& buffers_);


	void ClearPairBuffer();
	void CompressPair(const FastqRecord& record_);


	// pair matching
	//
	struct PairMatch
	{
		const FastqRecord* rec = nullptr;
		std::array<uint16, 4> sigPos;

		void Clear()
		{
			rec = nullptr;
			std::fill(sigPos.begin(), sigPos.end(), 0);
		}
	};

	struct LzPairMatch : public PairMatch
	{
		std::array<uint32, 4> signatures;

		void Clear()
		{
			PairMatch::Clear();
			std::fill(signatures.begin(), signatures.end(), 0);
		}
	};


	// pair encoding context
	//
	typedef TAdvancedContextCoder<8, 4> FlagContextCoder;

	struct DnaEncodersPE : public DnaEncoders
	{
		typedef TEncoder<FlagContextCoder> FlagEncoder;

		FlagEncoder* flagCoder;

		DnaEncodersPE()
			:	flagCoder(NULL)
		{}
	};

	DnaEncodersPE pairCtx;
	ReadMatchType curReadMatchTypeSe;

	std::deque<LzPairMatch*> pairHistory;
	//std::map<uint32, std::vector<const LzPairMatch*>> pairSignatureMap;
	std::unordered_map<uint32, std::vector<const LzPairMatch*>> pairSignatureMap;

	FastqCategorizerBase categorizer;

	// for dry run PE
	Buffer* dryFastqBuffer_2;
	BitMemoryWriter* dryFastqWriter_2;

};


class LzDecompressorPE : public LzDecompressorSE
{
public:
	using LzDecompressorSE::LzDecompressorSE;

	void Decompress(CompressedFastqBlock& compBin_,
					  std::vector<FastqRecord>& reads_,
					  FastqCompressedBin& fastqWorkBin_,
					  FastqChunk& dnaBuffer_);

private:
	enum ReadFlagsPE
	{
		ReadDifficultPE = 0,
		ReadIdenticalPE,
		ReadFullEncodePE,
		ReadFullExpensivePE,
		ReadShiftOnlyPE,
	};

	void DecompressMetaAndLinkRecords(std::vector<FastqRecord>& reads_, FastqChunk& chunk_);

	void SetupBufferMask(std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_);
	void StartDecoding(std::vector<DataChunk*>& buffers_);
	void EndDecoding();

	void DecompressRecord(FastqRecord& rec_, ReadFlags flag_,
						  const FastqRecord* auxRec_ = NULL, bool aux_ = false);

	void DecompressPair(FastqRecord& record_);


	// pair decoding context
	//
	typedef TAdvancedContextCoder<8, 4> FlagContextCoder;

	struct DnaDecodersPE : public DnaDecoders
	{
		typedef TDecoder<FlagContextCoder> FlagDecoder;

		FlagDecoder* flagCoder;

		DnaDecodersPE()
			:	flagCoder(NULL)
		{}
	};

	ReadMatchType curReadMatchTypeSe;

	DnaDecodersPE pairCtx;
	std::deque<const FastqRecord*> pairHistory;
};



class RawCompressorPE : public RawCompressorSE
{
public:
	RawCompressorPE(const CompressorParams& params_,
					const QualityCompressionData& globalQuaData_,
					const FastqRawBlockStats::HeaderStats& headData_,
					const CompressorAuxParams& auxParams_);
	~RawCompressorPE();


private:
	void DryRunCompress(const FastqRecord& record_, BitMemoryWriter* dryFastqWriter_);

	void CompressReadSequence(const FastqRecord& rec_);
	void CompressReadQuality(const FastqRecord& rec_);

	void StartEncoding(std::vector<DataChunk*>& buffers_);
	void EndEncoding(std::vector<DataChunk*>& buffers_);


	// for dry run PE
	Buffer* dryFastqBuffer_2;
	BitMemoryWriter* dryFastqWriter_2;
};


class RawDecompressorPE : public RawDecompressorSE
{
public:
	using RawDecompressorSE::RawDecompressorSE;

private:
	void DecompressRecords(std::vector<FastqRecord>& reads_, FastqChunk& Buffer_);

	void DecompressReadSequence(FastqRecord& rec_);
	void DecompressReadQuality(FastqRecord& rec_);

};




/**
 * FASTQ reads compressor/decompressor proxies
 *
 */
class FastqCompressor
{
public:
	FastqCompressor(const CompressorParams& params_,
					const QualityCompressionData& globalQuaData_,
					const FastqRawBlockStats::HeaderStats& headData_,
					const CompressorAuxParams& auxParams_ = CompressorAuxParams());
	~FastqCompressor();

	void Compress(const std::vector<FastqRecord>& reads_,
					 PackContext& packCtx_,
					 uint32 minimizerId_,
					 uint64 rawDnaStreamSize_,
					 FastqCompressedBin& fastqWorkBin_,
					 CompressedFastqBlock &compBin_);


	uint64 DecompressNext(CompressedFastqBlock& compBin_,
					   std::vector<FastqRecord>& reads_,
					   FastqCompressedBin& fastqWorkBin_,
					   FastqChunk& dnaBuffer_,
					   uint64 recStartIdx_ = 0);
private:
	const CompressorParams& params;
	const QualityCompressionData& globalQuaData;
	const FastqRawBlockStats::HeaderStats& headData;
	const CompressorAuxParams& auxParams;

	// TODO: refactor as it will have 2 sets of PPMd buffers  -- USE a pointer to PPMd compressor instead of object inside
	//
	LzCompressorSE* lzCompressor;
	RawCompressorSE* rawCompressor;
};


class FastqDecompressor
{
public:
	FastqDecompressor(const CompressorParams& params_,
					const QualityCompressionData& globalQuaData_,
					const FastqRawBlockStats::HeaderStats& headData_,
					const CompressorAuxParams& auxParams_ = CompressorAuxParams());
	~FastqDecompressor();

	// old routines
	//
	void Decompress(CompressedFastqBlock& compBin_,
					   std::vector<FastqRecord>& reads_,
					   FastqCompressedBin& fastqWorkBin_,
					   FastqChunk& dnaBuffer_);



private:
	const CompressorParams& params;
	const QualityCompressionData& globalQuaData;
	const FastqRawBlockStats::HeaderStats& headData;
	const CompressorAuxParams& auxParams;

	LzDecompressorSE* lzDecompressor;
	RawDecompressorSE* rawDecompressor;
};


#if 0
class FastqDecompressorDyn
{
public:
	FastqDecompressorDyn(const CompressorParams& params_,
					const QualityCompressionData& globalQuaData_,
					const FastqRawBlockStats::HeaderStats& headData_,
					const CompressorAuxParams& auxParams_ = CompressorAuxParams());
	~FastqDecompressorDyn();


	// new routines
	//
	void StartDecompress(CompressedFastqBlock& compBin_,
					   FastqCompressedBin& fastqWorkBin_);

	uint64 DecompressNext(std::vector<FastqRecord>& reads_,
						  FastqChunk& dnaBuffer_,
						  uint64 recStartIdx_ = 0);

	void FinishDecompress();



private:
	const CompressorParams& params;
	const QualityCompressionData& globalQuaData;
	const FastqRawBlockStats::HeaderStats& headData;
	const CompressorAuxParams& auxParams;

	LzDecompressorSE* lzDecompressor;
	RawDecompressorSE* rawDecompressor;
};
#endif


#endif // H_BINCOMPRESSOR
