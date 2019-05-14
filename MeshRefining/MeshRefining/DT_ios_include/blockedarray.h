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
		wordsPerBlock = 0;		/* ÿһ����ڴ��С����λ���ֽ� */
		wordsPerItem = 0;		/* ÿһ��Ԫ�ص��ڴ��С����λ���ֽ� */
		itemCntPerBlock = 0;	/* ÿһ��������ٸ�Ԫ�أ�2��ָ�� */
	
		activeItemCnt = 0;		/* ��ռ���˶��ٸ�Ԫ�صĿռ� */
		allocatedItemCnt = 0;	/* �ѷ����ڴ湻����Ԫ��ʹ�� */
	
		activeBlockCnt = 0;		/* ��ʹ�õĿ����Ŀ */
		allocatedBlockCnt = 0;	/* �ѷ����ڴ�Ŀ����Ŀ */	 
		allowedBlockCnt = 0;	/* ��������Ŀ����Ŀ */ 
	
		blocks = NULL;			/* ���п��ָ�� */
#if 0
		activeblock = NULL;		/* ����ʹ�õ���Ч�� */
		firstblock = NULL;		/* ��һ���ָ�� */
		lastvalidblock = NULL;	/* ���һ����Ч��(�ѷ������ڴ�����һ��) */
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

		activeItemCnt = 0;		/* ��ռ���˶��ٸ�Ԫ�صĿռ� */
		allocatedItemCnt = 0;	/* �ѷ����ڴ湻����Ԫ��ʹ�� */
	
		activeBlockCnt = 0;					/* ��ʹ�õĿ����Ŀ */
		allocatedBlockCnt = 0;				/* �ѷ����ڴ�Ŀ����Ŀ */	 
		allowedBlockCnt = INIT_BLOCK_CNT;	/* ��������Ŀ����Ŀ */ 

#if 0	
		activeblock = NULL;		/* ����ʹ�õ���Ч�� */
		firstblock = NULL;		/* ��һ���ָ�� */
		lastvalidblock = NULL;	/* ���һ����Ч��(�ѷ������ڴ�����һ��) */
		lastblock = allowedBlockCnt <= 0 ? NULL : blocks[allowedBlockCnt - 1];
#endif

		return 0;	/* �ɹ� */
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
		{/* ��Ҫ�����µ��ڴ� */
			newBlockCnt = (nTotalItems - allocatedItemCnt) / itemCntPerBlock + 
				((nTotalItems - allocatedItemCnt) % itemCntPerBlock > 0 ? 1 : 0);
			
			/* �����ǰ���õĿ鲻�������Է����µ��ڴ� */
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

			/* Ϊ�洢�µ�Ԫ�أ������µĿ� */
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
		
		wordsPerBlock = 0;		/* ÿһ����ڴ��С����λ���ֽ� */
		wordsPerItem = 0;		/* ÿһ��Ԫ�ص��ڴ��С����λ���ֽ� */
		itemCntPerBlock = 0;	/* ÿһ��������ٸ�Ԫ�أ�2��ָ�� */
	
		activeItemCnt = 0;		/* ��ռ���˶��ٸ�Ԫ�صĿռ� */
		allocatedItemCnt = 0;	/* �ѷ����ڴ湻����Ԫ��ʹ�� */
	
		activeBlockCnt = 0;		/* ��ʹ�õĿ����Ŀ */
		allocatedBlockCnt = 0;	/* �ѷ����ڴ�Ŀ����Ŀ */	 
		allowedBlockCnt = 0;	/* ��������Ŀ����Ŀ */ 
	
		blocks = NULL;			/* ���п��ָ�� */
	}

	/* ��ǰ�ƽ���ռ��һ���µĵ�Ԫ��λ */
	int occupyTailElem()
	{
		/* ��ǰ�ƽ�һ���µ�Ԫ */
		if (activeItemCnt + 1 > allocatedItemCnt)
			return -1; /* ��Ҫ�����µ��ڴ�ſ��� */
		else
		{
			activeItemCnt++; /* �µĵ�Ԫ */
			/* ����µ�Ԫ��Խ����һ�����У�����Ҫ���������Ϣ */
			if ((activeItemCnt - 1) % itemCntPerBlock == 0)
			{
				activeBlockCnt++;
#if 0
				activeblock = activeBlockCnt <= 0 ? NULL : blocks[activeBlockCnt - 1];
#endif
			}
			return activeItemCnt - 1; /* ��ǰʹ�õ�Item���� */
		}
	}

	int leaveTailElem()
	{
		assert(activeItemCnt <= allocatedItemCnt);
		assert(activeItemCnt > 0);	

		activeItemCnt--; 
		/* ����µ�Ԫ��Խ����һ�����У�����Ҫ���������Ϣ */
		if ((activeItemCnt + 1) % itemCntPerBlock == 0)
		{
			activeBlockCnt--;
#if 0
			activeblock = activeBlockCnt <= 0 ? NULL : blocks[activeBlockCnt - 1];
#endif
		}

		return activeItemCnt - 1; /* ��ǰʹ�õ�Item���� */
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

	/* ������������Ԫ��ֵ�����ز����� */
	Type& operator[](const int idx) const
	{
		/* ��Ч������ֵ��[0, activeItemCnt)֮�䣬Ϊ�˿��ٷ��ʣ�������� */
		return blocks[(idx >> log2_itemCntPBlk)][(idx) & itemCntPerBlockMinusOne]; 
	}

protected:
	int wordsPerBlock;			/* ÿһ����ڴ��С����λ���ֽ� */
	int wordsPerItem;			/* ÿһ��Ԫ�ص��ڴ��С����λ���ֽ� */
	int itemCntPerBlock;		/* ÿһ��������ٸ�Ԫ�� */
	int itemCntPerBlockMinusOne;/* itemCntPerBlock - 1������ͨ��idx���ٷ���Ԫ�� */
	int log2_itemCntPBlk;		/* ÿһ��������ٸ�Ԫ�أ�2��ָ��ֵ */

	int activeItemCnt;		/* ��ռ���˶��ٸ�Ԫ�صĿռ� */
	int allocatedItemCnt;	/* �ѷ����ڴ湻����Ԫ��ʹ�� */
	
	int activeBlockCnt;		/* ��ʹ�õĿ����Ŀ */
	int allocatedBlockCnt;	/* �ѷ����ڴ�Ŀ����Ŀ */	 
	int allowedBlockCnt;	/* ��������Ŀ����Ŀ */ 

	Type **blocks;			/* ���п��ָ�� */
#if 0
	Type *activeblock;		/* ����ʹ�õ���Ч�� */
	Type *firstblock;		/* ��һ���ָ�� */
	Type *lastvalidblock;	/* ���һ����Ч��(�ѷ������ڴ�����һ��) */
	Type *lastblock;		/* ���һ�� (�����ǿ�ָ��) */
#endif
};

#undef Type
#undef BlockedArray_Type