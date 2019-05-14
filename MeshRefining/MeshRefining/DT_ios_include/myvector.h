#ifndef __my_vector_h__
#define __my_vector_h__

#include <assert.h>

#define VECTOR_SIZE_ADD_RATIO 1.2
#define VECTOR_SIZE_MIN_NUM 8192 /* 2**13 */

template <class T>  
class MyVector{  
public:  
      
    typedef T * iterator;  
    typedef T * pointer;  
    typedef T & reference;  
      
    MyVector():start(0),finish(0),end_of_storage(0){}  
    MyVector(int n,const reference val)//申请n单位空间，并用val初始化  
    {  
        start = (T*) malloc(sizeof(T)*n);
		if (!start)
		{
			printf("Not enough memory.\n");
			exit(1);
		}
        finish = start;  
        end_of_storage = start + n;  
        while(n--) { *finish++=val; }  
    }  
    MyVector(int n)//申请n单位空间，初始化为0  
    {  
		start = (T*) malloc(sizeof(T)*n);
		if (!start)
		{
			printf("Not enough memory.\n");
			exit(1);
		} 
        finish = start;  
        end_of_storage = start + n;  
        while(n--) { *finish++=0; }  
    }  
      
    ~MyVector(){  
        iterator i;  
        for(i=start;i<end_of_storage;i++) i->~T();  
    }  
      
    iterator begin() const { return start; }  
    iterator end() const { return finish; }  
      
    //注意size的定义，返回值是end - begin，也就是说，finish指向的是最后一个元素后面的一个地址  
    int size() const { return int(finish-start); }  
      
    //void resize(int new_size,const reference x);//重新定义空间大小，而且完成初始化，finish指针在空间尾端  
    //void resize(int new_size);  
    //void reserve(int new_size);//重新定义空间大小，但是不改变finish指针  
    int capacity() const { return end_of_storage-start; }  
      
    bool empty() const { return start == finish; }  
      
    reference operator[](int n) { return *(start+n); }  
    reference front() { return *start; }  
    reference back() { return *(finish-1); }  
      
      
    void push_back(const reference x){
		int maxSize, newSize;
        if(finish == end_of_storage)//未满载  
        {  //满载，需要重新分配空间  
			maxSize = capacity();
			newSize = maxSize*VECTOR_SIZE_ADD_RATIO > VECTOR_SIZE_MIN_NUM ? 
				maxSize*VECTOR_SIZE_ADD_RATIO : VECTOR_SIZE_MIN_NUM;
			assert(newSize > maxSize + 1);
			start = (T*)realloc(start, sizeof(T)*newSize);
			if (!start)
			{
				printf("Not enough memory.\n");
				exit(1);
			}
			finish = start + maxSize;
			end_of_storage = start + newSize;
		}
		*finish++ = x;  
     //   finish++; 
    }  
      
    void pop_back()//删除尾端元素  
    {  
        finish--;  
        finish->~T();  
    }  
      
    void erase(iterator first,iterator last)//清除[first,last)的所有元素  
    {  
        int j = end()-last;  
        for(int i=0;i<j;i++){  
            *(first+i) = *(last+i);  
        }  
        while(end()>first+j){  
            pop_back();  
        }  
    }  
      
    void erase(iterator position)//清除指定位置元素  
    {  
        erase(position, position+1);  
    }  
      
    void clear(){//清空vector  
 //       erase(begin(),end());  
		finish = start;
    }  
      
    void insert(iterator position,int n,const reference x);//从position开始插入n个元素，每个元素的初始值为x  
    //void insert_aux(iterator position, const reference x);//重新分配空间，并完成已有成员的搬移  
      
private:  
    iterator start;  
    iterator finish;  
    iterator end_of_storage;  
};  
  
#if 0
template <class T>  
void MyVector<T>::insert(iterator position,int n,const reference x)//从position开始插入n个元素，每个元素的初始值为x  
{  
    /* 
    STL的vector<class type>::insert(type * , int , const type &) 作用：在指定位置插入元素，有可能引起扩容 
    当插入导致finish 指针大于 end_ofstorage 时，就需要扩容 
    */  
      
    iterator i = start;//old_start  
    iterator new_finish;//不扩容的话指向old_start，扩容的话指向new_start  
    iterator old_start = start,old_finish = finish;  
    bool needToDestory = false;  
    if(finish+n > end_of_storage){  
  
        needToDestory = true;  
  
        const int old_size = size();  
        const int len = old_size + n;  
  
        start = new T[len];  
        new_finish = start;  
        end_of_storage = start + len;  
    }  
    else {  
        new_finish = start;  
        end_of_storage = finish + n;  
    }  
      
    /*原始数据的搬移*/  
    while(i<position) { *new_finish++ = *i++; }  
    /*插入部分的初始化*/  
    while(n--) { *new_finish++ = x; }  
    /*原始数据的搬移*/  
    while(i<finish) { *new_finish++ = *i++; }  
    finish = new_finish;  
  
    /*这里需要释放原有空间*/  
    if(needToDestory==true){  
        while(old_start != old_finish){  
            old_start->~T();  
            old_start++;  
        }  
    }  
}  
#endif

#endif /* #ifndef __my_vector_h__ */  
