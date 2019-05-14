#include <memory.h>
#include <assert.h>

/* ---------------------------------------------------------------------------------------------
 * blocks -->-----------------------<--------firstblock
 *          |   0                  |
 *          -----------------------
 *          |   1                  |
 *          -----------------------
 *          |   ...                |
 *          ------------------------<--------activeblock
 *          | activeBlockCnt  - 1    |
 *          ----------------------- 
 *          |   ...                |
 *          ------------------------
 *          | allocatedBlockCnt - 1|
 *          -----------------------<---------lastvalidblock
 
 *          |   ...                |
 *          ------------------------
 *          | allowedBlockCnt - 1  |
 *          -----------------------<---------lastblock
 * ------------------------------------------------------------------------------------------*/
#define INIT_BLOCK_CNT	1024
#define BLOCK_ADD_CNT	128

class BlockedArray_Type
{
public:
	BlockedArray_Type() 
	{
		wordsPerBlock = 0;		/* 每一块的内存大小，单位：字节 */
		wordsPerItem = 0;		/* 每一个元素的内存大小，单位：字节 */
		itemCntPerBlock = 0;	/* 每一块包含多少个元素，2的指数 */
	
		activeItemCnt = 0;		/* 已占用了多少个元素的空间 */
		allocatedItemCnt = 0;	/* 已分配内存够多少元素使用 */
	
		activeBlockCnt = 0;		/* 已使用的块的数目 */
		allocatedBlockCnt = 0;	/* 已分配内存的块的数目 */	 
		allowedBlockCnt = 0;	/* 允许的最多的块的数目 */ 
	
		blocks = NULL;			/* 所有块的指针 */
#if 0
		activeblock = NULL;		/* 正在使用的有效块 */
		firstblock = NULL;		/* 第一块的指针 */
		lastvalidblock = NULL;	/* 最后一个有效块(已分配有内存的最后一块) */
#endif
	}
	~BlockedArray_Type() 
	{
		freeMemory();
	}

	int init(int log2ItemPBlk) 
	{
		log2_itemCntPBlk = log2ItemPBlk;
		itemCntPerBlock = 1 << log2_itemCntPBlk;	
		itemCntPerBlockMinusOne = itemCntPerBlock - 1;
		wordsPerItem = sizeof(Type);
		wordsPerBlock = wordsPerItem * itemCntPerBlock;

		blocks = (Type**) malloc(sizeof(Type*)*INIT_BLOCK_CNT);
		if (!blocks)
		{
			printf("Not enough memory.\n");
			exit(1);
		}
		memset(blocks, 0, sizeof(Type*)*INIT_BLOCK_CNT);

		activeItemCnt = 0;		/* 已占用了多少个元素的空间 */
		allocatedItemCnt = 0;	/* 已分配内存够多少元素使用 */
	
		activeBlockCnt = 0;					/* 已使用的块的数目 */
		allocatedBlockCnt = 0;				/* 已分配内存的块的数目 */	 
		allowedBlockCnt = INIT_BLOCK_CNT;	/* 允许的最多的块的数目 */ 

#if 0	
		activeblock = NULL;		/* 正在使用的有效块 */
		firstblock = NULL;		/* 第一块的指针 */
		lastvalidblock = NULL;	/* 最后一个有效块(已分配有内存的最后一块) */
		lastblock = allowedBlockCnt <= 0 ? NULL : blocks[allowedBlockCnt - 1];
#endif

		return 0;	/* 成功 */
	}

	int clean()
	{
		int i;
		for (i = 0; i < allocatedBlockCnt; i++)
			memset(blocks[i], 0, wordsPerBlock);
		return 1;
	}

	int allocMemory(int nTotalItems)
	{
		int i, newBlockCnt = 0, newAllowedBlockCnt = allowedBlockCnt;
		Type **blocks_new = NULL;

		if (nTotalItems > allocatedItemCnt)
		{/* 需要分配新的内存 */
			newBlockCnt = (nTotalItems - allocatedItemCnt) / itemCntPerBlock + 
				((nTotalItems - allocatedItemCnt) % itemCntPerBlock > 0 ? 1 : 0);
			
			/* 如果当前可用的块不够，则尝试分配新的内存 */
			while (newAllowedBlockCnt < newBlockCnt + allocatedBlockCnt)
				newAllowedBlockCnt += BLOCK_ADD_CNT;
			if (newAllowedBlockCnt > allowedBlockCnt)
			{
				blocks_new = (Type**) realloc(blocks, sizeof(Type*)*newAllowedBlockCnt);
				if (!blocks_new)
				{
					printf("Not enough memory.\n");
					exit(1);
				}
				memset(blocks_new + allowedBlockCnt, 0, sizeof(Type*)*(newAllowedBlockCnt - allowedBlockCnt));
				if (blocks != blocks_new)
				{
					blocks = blocks_new;
#if 0
					activeblock = activeBlockCnt <= 0 ? NULL : blocks[activeBlockCnt - 1];
					firstblock = blocks[0];
					lastvalidblock = allocatedBlockCnt <= 0 ? NULL : blocks[allocatedBlockCnt - 1];
#endif
					
				}

				allowedBlockCnt = newAllowedBlockCnt;
#if 0
				lastblock = allowedBlockCnt <= 0 ? NULL : blocks_new[allowedBlockCnt - 1];
#endif
			}

			/* 为存储新的元素，分配新的块 */
			for (i = 0; i < newBlockCnt; i++)
			{
				blocks[i+allocatedBlockCnt] = (Type*) malloc(wordsPerBlock);  
				if (!blocks[i+allocatedBlockCnt])
				{
					printf("Not enough memory.\n");
					exit(1);
				}
				memset(blocks[i+allocatedBlockCnt], 0, wordsPerBlock);
			}
			allocatedBlockCnt += newBlockCnt;
			allocatedItemCnt += newBlockCnt * itemCntPerBlock;

#if 0
			lastvalidblock = allocatedBlockCnt <= 0 ? NULL : blocks[allocatedBlockCnt - 1];
#endif
		}

		return allocatedItemCnt;
	}

	void freeMemory()
	{
		int i;
		for (i = 0; i < allocatedBlockCnt; i++)
			if (blocks[i])
				free(blocks[i]);
		if (blocks)
			free(blocks);
		
		wordsPerBlock = 0;		/* 每一块的内存大小，单位：字节 */
		wordsPerItem = 0;		/* 每一个元素的内存大小，单位：字节 */
		itemCntPerBlock = 0;	/* 每一块包含多少个元素，2的指数 */
	
		activeItemCnt = 0;		/* 已占用了多少个元素的空间 */
		allocatedItemCnt = 0;	/* 已分配内存够多少元素使用 */
	
		activeBlockCnt = 0;		/* 已使用的块的数目 */
		allocatedBlockCnt = 0;	/* 已分配内存的块的数目 */	 
		allowedBlockCnt = 0;	/* 允许的最多的块的数目 */ 
	
		blocks = NULL;			/* 所有块的指针 */
	}

	/* 往前推进，占有一个新的单元空位 */
	int occupyTailElem()
	{
		/* 往前推进一个新单元 */
		if (activeItemCnt + 1 > allocatedItemCnt)
			return -1; /* 需要分配新的内存才可以 */
		else
		{
			activeItemCnt++; /* 新的单元 */
			/* 如果新单元跨越到下一个块中，则需要更新相关信息 */
			if ((activeItemCnt - 1) % itemCntPerBlock == 0)
			{
				activeBlockCnt++;
#if 0
				activeblock = activeBlockCnt <= 0 ? NULL : blocks[activeBlockCnt - 1];
#endif
			}
			return activeItemCnt - 1; /* 当前使用的Item数量 */
		}
	}

	int leaveTailElem()
	{
		assert(activeItemCnt <= allocatedItemCnt);
		assert(activeItemCnt > 0);	

		activeItemCnt--; 
		/* 如果新单元跨越到下一个块中，则需要更新相关信息 */
		if ((activeItemCnt + 1) % itemCntPerBlock == 0)
		{
			activeBlockCnt--;
#if 0
			activeblock = activeBlockCnt <= 0 ? NULL : blocks[activeBlockCnt - 1];
#endif
		}

		return activeItemCnt - 1; /* 当前使用的Item数量 */
	}

	int getArraySize()
	{
		return activeItemCnt;
	}

	void setArraySize(int arraySize)
	{
		assert(arraySize < allocatedItemCnt);
		activeItemCnt = arraySize;
		activeBlockCnt = activeItemCnt / itemCntPerBlock + (activeItemCnt % itemCntPerBlock > 0 ? 1 : 0);
	}

	int getArrayCapacity()
	{
		return allocatedItemCnt;
	}

	/* 根据索引访问元素值，重载操作符 */
	Type& operator[](const int idx) const
	{
		/* 有效的索引值在[0, activeItemCnt)之间，为了快速访问，不做检查 */
		return blocks[(idx >> log2_itemCntPBlk)][(idx) & itemCntPerBlockMinusOne]; 
	}

protected:
	int wordsPerBlock;			/* 每一块的内存大小，单位：字节 */
	int wordsPerItem;			/* 每一个元素的内存大小，单位：字节 */
	int itemCntPerBlock;		/* 每一块包含多少个元素 */
	int itemCntPerBlockMinusOne;/* itemCntPerBlock - 1，用于通过idx快速访问元素 */
	int log2_itemCntPBlk;		/* 每一块包含多少个元素，2的指数值 */

	int activeItemCnt;		/* 已占用了多少个元素的空间 */
	int allocatedItemCnt;	/* 已分配内存够多少元素使用 */
	
	int activeBlockCnt;		/* 已使用的块的数目 */
	int allocatedBlockCnt;	/* 已分配内存的块的数目 */	 
	int allowedBlockCnt;	/* 允许的最多的块的数目 */ 

	Type **blocks;			/* 所有块的指针 */
#if 0
	Type *activeblock;		/* 正在使用的有效块 */
	Type *firstblock;		/* 第一块的指针 */
	Type *lastvalidblock;	/* 最后一个有效块(已分配有内存的最后一块) */
	Type *lastblock;		/* 最后一块 (可能是空指针) */
#endif
};

#undef Type
#undef BlockedArray_Type