/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2007��10��11��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 11, 2007
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953165
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��gridinterface.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ�����������������û���
 *     Define a set of data reference structures for grids
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
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
