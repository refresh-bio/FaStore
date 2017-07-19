/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../fastore_bin/Globals.h"

#include <vector>
#include <map>

#include "RebinOperator.h"

#include "../fastore_bin/BinFile.h"
#include "../fastore_bin/FastqCategorizer.h"
#include "NodesPacker.h"
#include "DnaRebalancer.h"


void BinBalancer::Run()
{
	const bool pairedEnd = binConfig.archiveType.readType == ArchiveType::READ_PE;

	int64 partId = 0;

	std::unique_ptr<IFastqNodesPackerDyn> packer(!pairedEnd
	  ? (IFastqNodesPackerDyn*)(new FastqNodesPackerDynSE(binConfig))
	  : (IFastqNodesPackerDyn*)(new FastqNodesPackerDynPE(binConfig)));

	DnaRebalancer rebalancer(binConfig.minimizer, balanceParams, pairedEnd);

	RebinWorkBuffer binBuffer;
	BinaryBinBlock* inPart = NULL;

	FastqRecordBinStats stats;
	IFastqChunkCollection tmpChunks;

	while (inPartsQueue->Pop(partId, inPart))
	{
		ASSERT(inPart->metaSize > 0);
		ASSERT(inPart->dnaSize > 0);
		ASSERT(inPart->rawDnaSize > 0);
		ASSERT(inPart->signature != 0);

		// only process the non-parity bins
		//
		const uint32 signatureId = inPart->signature;

		BinaryBinBlock* outPart = NULL;


		if (BinBalanceParameters::IsSignatureValid(signatureId, balanceParams.signatureParity))
		{
			packer->UnpackFromBin(*inPart,
								 binBuffer.reads,
								 *binBuffer.rebinCtx.graph,
								 stats,
								 tmpChunks,
								 false);

			ASSERT(inPart->rawDnaSize > 0);

			const uint64 inRawReadsCount = binBuffer.reads.size();

			rebalancer.Rebalance(binBuffer.rebinCtx,
								 binBuffer.nodesMap,
								 inPart->signature);

			// reclaim the used memory from the part
			//
			inPart->Clear();

			inPartsPool->Release(inPart);
			inPart = NULL;

			outPartsPool->Acquire(outPart);
			packer->PackToBins(binBuffer.nodesMap, *outPart);

			uint64 outRawReadsCount = 0;
			for (const auto& desc : outPart->descriptors)
				outRawReadsCount += desc.second.recordsCount;
			ASSERT(outRawReadsCount == inRawReadsCount);
		}
		else
		{
			// TODO: optimize merger: do not unpack, to save memory: just merge raw data or
			// change header of the block to being a multi-part bin
			//
			if (inPart->auxDescriptors.size() > 1)
			{
				// such case should not happen often --> only when re-binning from stage /n to /n+e, e > 1

				packer->UnpackFromBin(*inPart,
									 binBuffer.reads,
									 *binBuffer.rebinCtx.graph,
									 stats,
									 tmpChunks,
									 false);

				const uint64 inRawReadsCount = binBuffer.reads.size();

				// reclaim the used memory from the part
				//
				inPart->Clear();

				inPartsPool->Release(inPart);
				inPart = NULL;

				outPartsPool->Acquire(outPart);
				packer->PackToBin(*binBuffer.rebinCtx.graph, *outPart, signatureId);

				uint64 outRawReadsCount = 0;
				for (const auto& desc : outPart->auxDescriptors)
					outRawReadsCount += desc.recordsCount;
				ASSERT(outRawReadsCount == inRawReadsCount);
			}
			else
			{
				// just swap the parts
				//
				outPartsPool->Acquire(outPart);
				inPart->Swap(*outPart);

				inPart->Clear();
				inPartsPool->Release(inPart);
				inPart = NULL;
			}
		}

		// reclaim used memory
		//
		binBuffer.Reset();

#if EXTRA_MEM_OPT
		tmpChunks.Clear();
#endif

		outPartsQueue->Push(partId, outPart);
	}

	outPartsQueue->SetCompleted();
}
