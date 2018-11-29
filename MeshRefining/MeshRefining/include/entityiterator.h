/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2007年10月11日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 11, 2007
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953165
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：gridinterface.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对网格的数据引用机制
 *     Define a set of data reference structures for grids
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __eemas_entityiterator_h__
#define __eemas_entityiterator_h__

#ifdef EEMAS_USE_OLD_STD_HEADERS
# include "vector.h"
# include <list.h>
#else
# include <vector>
# include <list>
#endif
# include <assert.h>

typedef void* EntityHandle;

typedef std::vector<EntityHandle> EntityHandleArray;
typedef std::list<EntityHandle> EntityHandleList;

class EntityIterator
{
public:
	virtual ~EntityIterator() {}
	
	//! Moves the iterator back to the first
	//! entity in the list.
    virtual void restart() = 0;
	
	//! *iterator.  Return the handle currently
	//! being pointed at by the iterator.
    virtual EntityHandle operator*() const = 0;
    
    //! ++iterator
    virtual void operator++() = 0;

	//! --iterator
    virtual void operator--() = 0;
    
    //! Returns false until the iterator has
    //! been advanced PAST the last entity.
    //! Once is_end() returns true, *iterator
    //! returns 0.
    virtual bool is_end() const = 0;
};

class ArrayIterator: public EntityIterator
{
private:
	EntityHandleArray *entity_array;
	unsigned int cur_index;

public:
    ArrayIterator(EntityHandleArray *e_array = nullptr)
	{
		entity_array = e_array;
		cur_index = 0;
	}

	virtual ~ArrayIterator() 
	{
	}
	
	//! Moves the iterator back to the first
	//! entity in the list.
    virtual void restart()
	{
		cur_index = 0;
	}

	virtual void to_end()
	{
		cur_index = entity_array->size()-1;
	}
	
	//! *iterator.  Return the handle currently
	//! being pointed at by the iterator.
    virtual EntityHandle operator*() const
	{
		assert(entity_array && entity_array->size() > cur_index);
		return (*entity_array)[cur_index];
	}
    
    //! ++iterator
    virtual void operator++()
	{
		assert(entity_array && entity_array->size() > cur_index);
		cur_index++;
	}

	//! --iterator
    virtual void operator--()
	{
		assert(entity_array);
		cur_index--;
	}
    
    //! Returns false until the iterator has
    //! been advanced PAST the last entity.
    //! Once IsEnd() returns true, *iterator
    //! returns 0.
    virtual bool is_end() const
	{
		assert(entity_array && entity_array->size() >= cur_index);
		return cur_index == entity_array->size();
	}

	virtual bool is_head() const
	{
		assert(entity_array && entity_array->size() >= (cur_index+1));
		return (cur_index + 1) == 0;
	}

	virtual int size() const
	{
		assert(entity_array);
		return entity_array->size();
	}
};

class ListIterator: public EntityIterator
{
private:
	EntityHandleList *entity_list;
	EntityHandleList::iterator it;

public:
    ListIterator(EntityHandleList *e_list = nullptr)
	{
		entity_list = e_list;
        if (e_list != nullptr)
			it = entity_list->begin();
	}

	virtual ~ListIterator() 
	{
	}
	
	//! Moves the iterator back to the first
	//! entity in the list.
    virtual void restart()
	{
        assert(nullptr != entity_list);
		it = entity_list->begin();
		assert(it != entity_list->end());
		
	}
	void reend()
	{
        assert(nullptr != entity_list);
		it = entity_list->end();
		it--;
	}
	void clear()
	{
        assert(nullptr != entity_list);
		it = entity_list->begin();
		entity_list->clear();
	
		
	}
	
	//! *iterator.  Return the handle currently
	//! being pointed at by the iterator.
    virtual EntityHandle operator*() const
	{
		assert(entity_list && it != entity_list->end());
		return *it;
	}
    
    //! ++iterator
    virtual void operator++()
	{
		assert(entity_list && it != entity_list->end());
		it++;
	}

	//! --iterator
    virtual void operator--()
	{
		assert(entity_list && it != entity_list->end());
		it--;
	}
    
    //! Returns false until the iterator has
    //! been advanced PAST the last entity.
    //! Once IsEnd() returns true, *iterator
    //! returns 0.
    virtual bool is_end() const
	{
		assert(entity_list);
		return it == entity_list->end();
	}
	bool is_begin() const
	{
		assert(entity_list);
		return it == entity_list->begin();
	}
};

#endif /* __eemas_entityiterator_h__ */
