/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef NODE_H
#define NODE_H

#include "Globals.h"
#include "FastqRecord.h"

#include <list>


/**
 * Represents matching information which can be attached
 * to a node
 */
struct NodeDecorator
{
	enum DecoratorGroupType
	{
		GROUP_NONE = 0,
		GROUP_EXACT_MATCHES,
		GROUP_CONTIG,
		GROUP_SUBTREE,
		GROUP_TRANSTREE			// only used when rebinning
	};

	uint8 type;					// instead of virtual function
	// 7B padding ...

	union DecoratorGroupData
	{
		ExactMatchesGroup* exactMatches;
		ContigDefinition* contig;
		GraphEncodingContext* subTree;
		TreeTransferDefinition* transTree;
	} group;


	NodeDecorator()
		:	type(GROUP_NONE)
	{
		group.exactMatches = NULL;
	}
};



/**
 * Data types used to represent the relation graph between
 * the matched FASTQ reads and to store the related information
 *
 */
struct ExactMatchesGroup
{
	std::vector<FastqRecord*> records;
};


struct ConsensusDefinition
{
	std::vector<char> sequence;
	std::vector<bool> variantPositions;
	std::pair<uint32, uint32> range;

	uint32 variantsCount;
	uint32 readLen;			// deprecated


	ConsensusDefinition()
		:	variantsCount(0)
		,	readLen(0)
	{}
};


struct ContigDefinition
{
	ConsensusDefinition consensus;
	std::vector<MatchNode*> nodes;
};


struct TreeTransferDefinition
{
	uint32 signatureId;
	uint32 recordsCount;
	uint16 mainSignaturePos;

	TreeTransferDefinition()
		:	signatureId(0)
		,	recordsCount(0)
		,	mainSignaturePos(0)
	{}
};



/**
 * Describes the connection graph between matched reads
 *
 *  TODO: shall we extinguish between the main graph encoding context
 *  and the auxillary ones when encoding subtrees?
 */
struct GraphEncodingContext
{
	// main context
	//
	uint32 signatureId;
	int32 mainSignaturePos;
	std::array<char, MAX_SIGNATURE_LEN> signature;		// only for DEBUG


	// those needs to be pointers as when reading, we are directly
	// linking them with nodes as decorators
	std::vector<MatchNode> nodes;
	std::vector<GraphEncodingContext*> subTrees;
	std::vector<ExactMatchesGroup*> exactMatches;

	// encoders
	//
	GraphEncodingContext()
		:	signatureId(0)
		,	mainSignaturePos(0)
	{}

	~GraphEncodingContext()
	{
		Clear();
	}

	GraphEncodingContext* CreateSubTreeGroup()
	{
		subTrees.push_back(new GraphEncodingContext());
		return subTrees.back();
	}

	ExactMatchesGroup* CreateExactMatchesGroup()
	{
		exactMatches.push_back(new ExactMatchesGroup());
		return exactMatches.back();
	}

	void Clear()
	{
		//TODO: shall we use shrink_to_fit() ???
		//
		signatureId = 0;
		mainSignaturePos = 0;

		nodes.clear();

		for (auto sg : subTrees)
			delete sg;
		subTrees.clear();

		for (auto eg : exactMatches)
			delete eg;
		exactMatches.clear();

#if EXTRA_MEM_OPT
		nodes.shrink_to_fit();
		subTrees.shrink_to_fit();
		exactMatches.shrink_to_fit();
#endif
	}
};


/**
 * Describes the connection graph between matched reads, but
 * in rebinning stage
 *
 */
struct RebinContext
{
	GraphEncodingContext* graph;
	const bool ownsGraph;

	std::vector<TreeTransferDefinition*> transTrees;
	std::vector<MatchNode*> rootNodes;


	RebinContext(GraphEncodingContext* graph_ = NULL)
		:	graph(graph_)
		,	ownsGraph(graph_ == NULL)
	{
		if (ownsGraph)
			graph = new GraphEncodingContext();
	}

	~RebinContext()
	{
		Clear();

		if (ownsGraph)
			delete graph;
	}


	TreeTransferDefinition* CreateTransTreeGroup()
	{
		transTrees.push_back(new TreeTransferDefinition());
		return transTrees.back();
	}

	void Clear()
	{
		for (auto sg : transTrees)
			delete sg;
		transTrees.clear();

		rootNodes.clear();

		if (ownsGraph)
			graph->Clear();

#if EXTRA_MEM_OPT
		transTrees.shrink_to_fit();
		rootNodes.shrink_to_fit();
#endif
	}
};


/**
 * Describes the connection graph between matched reads, but
 * in packing stage
 *
 */
struct PackContext
{
	GraphEncodingContext* graph;
	const bool ownsGraph;

	std::vector<ContigDefinition*> contigs;
	std::vector<MatchNode*> rootNodes;
	FastqRecordBinStats stats;						// TODO: we can add more useful stuff here

	PackContext(GraphEncodingContext* graph_ = NULL)
		:	graph(graph_)
		,	ownsGraph(graph_ == NULL)
	{
		if (ownsGraph)
			graph = new GraphEncodingContext();
	}

	~PackContext()
	{
		Clear();

		if (ownsGraph)
			delete graph;
	}

	ContigDefinition* CreateContigGroup()
	{
		contigs.push_back(new ContigDefinition());
		return contigs.back();
	}

	void Clear(bool onlyContigs_ = false)
	{
		for (auto cg : contigs)
			delete cg;
		contigs.clear();

		if (onlyContigs_)
			return;

		rootNodes.clear();

		stats.Clear();

		if (ownsGraph)
			graph->Clear();

#if EXTRA_MEM_OPT
		contigs.shrink_to_fit();
		rootNodes.shrink_to_fit();
#endif
	}
};


/**
 * A node used to represent matching information
 *
 */
struct MatchNode
{
	enum NodeType
	{
		TYPE_NONE = 0,
		TYPE_HARD,
		TYPE_LZ,
		TYPE_CONTIG_READ
	};

	enum NodeDecoratorFlags
	{
		FLAG_HAS_EXACT_MATCHES	= BIT(0),
		FLAG_ENCODES_CONTIG		= BIT(1),
		FLAG_ENCODES_SUBTREE	= BIT(2),
		FLAG_ENCODES_TRANSTREE	= BIT(3)
	};

	enum NodeMatchingFlags
	{
		MATCH_SHIFT_ONLY		= BIT(4)
	};

	//static const uint32 MaxBranchSize = (uint16)-1;


	// 8B
	uint8 type;
	uint8 flags;
	int16 shiftValue;
	int16 encodeCost;
	//int16 branchSize;

	// 24B
	// TODO: optimize
	FastqRecord* record;
	FastqRecord* lzRecord;
	MatchNode* parentNode;

	// 16B
        std::list<MatchNode*>* children;			// TODO: try to move here to vector
	std::vector<NodeDecorator>* decorators;		// TODO: merge decorators with children??

	friend void swap(MatchNode& lmn_, MatchNode& rmn_) noexcept;

	MatchNode()
		:	type(TYPE_NONE)
		,	flags(0)
		,	shiftValue(0)
		,	encodeCost(0)
		//,	branchSize(0)
		,	record(NULL)
		,	lzRecord(NULL)
		,	parentNode(NULL)
		,	children(NULL)
		,	decorators(NULL)
	{}

	MatchNode(MatchNode&& mn_)
		:	MatchNode()
	{
		swap(*this, mn_);
	}

	MatchNode& operator=(MatchNode&& mn_)
	{
		swap(*this, mn_);
		return *this;
	}

	// MatchNode is non-copyable class
	MatchNode(const MatchNode& mn_) = delete;
	MatchNode& operator=(const MatchNode&) = delete;


	~MatchNode()
	{
		Clear();
	}

	void Clear()
	{
		type = TYPE_NONE;
		flags = 0;
		shiftValue = 0;
		encodeCost = 0;

		record = NULL;
		lzRecord = NULL;
		parentNode = NULL;

		if (children != NULL)
			delete children;
		children = NULL;

		if (decorators != NULL)
			delete decorators;
		decorators = NULL;
	}

	// returns the tree size
	// WARN: O(log n) complexity !
	//
	uint64 Size() const
	{
		uint64 size = 1;
		if (HasChildren())
		{
			for (auto c : *children)
				size += c->Size();
		}
		return size;
	}


	// family management
	//
	void AddChild(MatchNode* child_)
	{
		ASSERT(child_->record != record);

		if (children == NULL)
                        children = new std::list<MatchNode*>();
		else
			ASSERT(std::find(children->begin(), children->end(), child_) == children->end());

		children->push_back(child_);
	}

	bool HasChildren() const
	{
		return children != NULL && !children->empty();
	}

	void RemoveChild(MatchNode* child_)
	{
		ASSERT(child_->record != record);
		ASSERT(HasChildren());

                auto ic = std::find(children->begin(), children->end(), child_);
                ASSERT(ic != children->end());
                children->erase(ic);

		if (children->empty())
		{
			delete children;
			children = NULL;
		}
	}

	void RemoveChildren()
	{
		//ASSERT(HasChildren());
		if (children != NULL)
			delete children;
		children = NULL;
	}


	// matching flags handling
	//
	bool IsSetFlag(uint32 flag_) const
	{
		return (flags & flag_) != 0;
	}

	void SetFlag(uint32 flag_, bool b_)
	{
		b_ ? (flags |= flag_) : (flags &= ~flag_);
	}

	bool HasNoMismatches() const
	{
		return IsSetFlag(MATCH_SHIFT_ONLY);
	}

	void SetNoMismatches(bool b_)
	{
		SetFlag(MATCH_SHIFT_ONLY, b_);
	}


	// decorators handling
	//
	NodeDecorator& AddDecorator()
	{
		if (decorators == NULL)
			decorators = new std::vector<NodeDecorator>();
		decorators->push_back(NodeDecorator());
		return decorators->back();
	}


	// exact matches group used after perfmirming LZ matching
	//
	bool HasExactMatches() const
	{
		return IsSetFlag(FLAG_HAS_EXACT_MATCHES);
	}

	void CreateExactMatchesGroup(ExactMatchesGroup* ems_)
	{
		ASSERT(!HasExactMatches());

		NodeDecorator& group = AddDecorator();
		group.type = NodeDecorator::GROUP_EXACT_MATCHES;
		group.group.exactMatches = ems_;

		SetFlag(FLAG_HAS_EXACT_MATCHES, true);
	}

	void RemoveExactMatches()
	{
		ASSERT(HasExactMatches());
		ASSERT(decorators != NULL);

		auto iems = begin(*decorators);
		while (iems != end(*decorators))
		{
			if (iems->type == NodeDecorator::GROUP_EXACT_MATCHES)
			{
				iems = decorators->erase(iems);
				break;
			}
			iems++;
		}
		ASSERT(iems != end(*decorators) || decorators->empty());

		if (decorators->empty())
		{
			delete decorators;
			decorators = NULL;
		}

		SetFlag(FLAG_HAS_EXACT_MATCHES, false);
	}

	ExactMatchesGroup* GetExactMatches() const
	{
		ASSERT(HasExactMatches());
		ASSERT(decorators != NULL);

		auto iems = begin(*decorators);
		while (iems != end(*decorators))
		{
			if (iems->type == NodeDecorator::GROUP_EXACT_MATCHES)
				break;
			iems++;
		}
		ASSERT(iems != end(*decorators));
		return iems->group.exactMatches;
	}

	void AddExactMatch(FastqRecord* rec_)
	{
		ASSERT(decorators != NULL);
		ASSERT(HasExactMatches());

		NodeDecorator* ems = NULL;
		for (auto& deco : *decorators)
		{
			if (deco.type == NodeDecorator::GROUP_EXACT_MATCHES)
			{
				ems = &deco;
				break;
			}
		}

		ASSERT(ems != NULL);
		ems->group.exactMatches->records.push_back(rec_);
	}


	// contig group used after semi-assembly step
	//
	bool HasContigGroup() const
	{
		return IsSetFlag(FLAG_ENCODES_CONTIG);
	}

	void AddContigGroup(ContigDefinition* contig_)
	{
		ASSERT(!HasContigGroup());

		NodeDecorator& group = AddDecorator();
		group.type = NodeDecorator::GROUP_CONTIG;
		group.group.contig = contig_;

		SetFlag(FLAG_ENCODES_CONTIG, true);
	}

	ContigDefinition* GetContigGroup() const
	{
		ASSERT(HasContigGroup());

		ContigDefinition* contig = NULL;
		for (auto& deco : *decorators)
		{
			if (deco.type == NodeDecorator::GROUP_CONTIG)
			{
				contig = deco.group.contig;
				break;
			}
		}

		ASSERT(contig != NULL);
		return contig;
	}


	// sub-tree group used after rebining
	//
	bool HasSubTreeGroup() const
	{
		return IsSetFlag(FLAG_ENCODES_SUBTREE);
	}

	void AddSubTreeGroup(GraphEncodingContext* tree_)
	{
		NodeDecorator& group = AddDecorator();
		group.type = NodeDecorator::GROUP_SUBTREE;
		group.group.subTree = tree_;

		SetFlag(FLAG_ENCODES_SUBTREE, true);
	}

	std::vector<GraphEncodingContext*> GetSubTrees() const
	{
		ASSERT(HasSubTreeGroup());

		std::vector<GraphEncodingContext*> trees;

		for (auto& deco : *decorators)
		{
			if (deco.type != NodeDecorator::GROUP_SUBTREE)
				continue;

			trees.push_back(deco.group.subTree);
		}

		ASSERT(trees.size() > 0);
		return trees;
	}

	void RemoveSubTrees()
	{
		ASSERT(HasSubTreeGroup());
		ASSERT(decorators != NULL);

		auto iems = begin(*decorators);
		while (iems != end(*decorators))
		{
			if (iems->type == NodeDecorator::GROUP_SUBTREE)
				iems = decorators->erase(iems);
			else
				iems++;
		}
		ASSERT(iems == end(*decorators));


		if (decorators->empty())
		{
			delete decorators;
			decorators = NULL;
		}

		SetFlag(FLAG_ENCODES_SUBTREE, false);
	}


	// transfer tree definition used when rebining
	//
	bool HasTransTreeGroup() const
	{
		return IsSetFlag(FLAG_ENCODES_TRANSTREE);
	}

	void AddTransTreeGroup(TreeTransferDefinition* tree_)
	{
		NodeDecorator& group = AddDecorator();
		group.type = NodeDecorator::GROUP_TRANSTREE;
		group.group.transTree = tree_;

		SetFlag(FLAG_ENCODES_TRANSTREE, true);
	}

	TreeTransferDefinition* GetTransTree() const
	{
		ASSERT(HasTransTreeGroup());

		TreeTransferDefinition* tree = NULL;
		for (auto& deco : *decorators)
		{
			if (deco.type == NodeDecorator::GROUP_TRANSTREE)
			{
				tree = deco.group.transTree;
				break;
			}
		}

		ASSERT(tree != NULL);
		return tree;
	}
};


inline void swap(MatchNode& lmn_, MatchNode& rmn_) noexcept
{
	std::swap(lmn_.type,		rmn_.type);
	std::swap(lmn_.flags,		rmn_.flags);
	std::swap(lmn_.shiftValue,	rmn_.shiftValue);
	std::swap(lmn_.encodeCost,	rmn_.encodeCost);
	std::swap(lmn_.record,		rmn_.record);
	std::swap(lmn_.lzRecord,	rmn_.lzRecord);
	std::swap(lmn_.parentNode,	rmn_.parentNode);
	std::swap(lmn_.children,	rmn_.children);
	std::swap(lmn_.decorators,	rmn_.decorators);
}


/**
 * Comparator used to sort the nodes (similarily as regular FASTQ reads)
 *
 */
template <>
struct TFastqComparator<const MatchNode*> : public IFastqComparator
{
	bool operator() (const MatchNode* n1_, const MatchNode* n2_)
	{
		return CompareReads(*n1_->record, *n2_->record);
	}
};

template <>
struct TFastqComparator<const MatchNode&> : public IFastqComparator
{
	bool operator() (const MatchNode& n1_, const MatchNode& n2_)
	{
		return CompareReads(*n1_.record, *n2_.record);
	}
};


#endif // NODE_H
