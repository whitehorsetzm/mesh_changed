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
    MyVector(int n,const reference val)//����n��λ�ռ䣬����val��ʼ��  
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
    MyVector(int n)//����n��λ�ռ䣬��ʼ��Ϊ0  
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
      
    //ע��size�Ķ��壬����ֵ��end - begin��Ҳ����˵��finishָ��������һ��Ԫ�غ����һ����ַ  
    int size() const { return int(finish-start); }  
      
    //void resize(int new_size,const reference x);//���¶���ռ��С��������ɳ�ʼ����finishָ���ڿռ�β��  
    //void resize(int new_size);  
    //void reserve(int new_size);//���¶���ռ��С�����ǲ��ı�finishָ��  
    int capacity() const { return end_of_storage-start; }  
      
    bool empty() const { return start == finish; }  
      
    reference operator[](int n) { return *(start+n); }  
    reference front() { return *start; }  
    reference back() { return *(finish-1); }  
      
      
    void push_back(const reference x){
		int maxSize, newSize;
        if(finish == end_of_storage)//δ����  
        {  //���أ���Ҫ���·���ռ�  
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
      
    void pop_back()//ɾ��β��Ԫ��  
    {  
        finish--;  
        finish->~T();  
    }  
      
    void erase(iterator first,iterator last)//���[first,last)������Ԫ��  
    {  
        int j = end()-last;  
        for(int i=0;i<j;i++){  
            *(first+i) = *(last+i);  
        }  
        while(end()>first+j){  
            pop_back();  
        }  
    }  
      
    void erase(iterator position)//���ָ��λ��Ԫ��  
    {  
        erase(position, position+1);  
    }  
      
    void clear(){//���vector  
 //       erase(begin(),end());  
		finish = start;
    }  
      
    void insert(iterator position,int n,const reference x);//��position��ʼ����n��Ԫ�أ�ÿ��Ԫ�صĳ�ʼֵΪx  
    //void insert_aux(iterator position, const reference x);//���·���ռ䣬��������г�Ա�İ���  
      
private:  
    iterator start;  
    iterator finish;  
    iterator end_of_storage;  
};  
  
#if 0
template <class T>  
void MyVector<T>::insert(iterator position,int n,const reference x)//��position��ʼ����n��Ԫ�أ�ÿ��Ԫ�صĳ�ʼֵΪx  
{  
    /* 
    STL��vector<class type>::insert(type * , int , const type &) ���ã���ָ��λ�ò���Ԫ�أ��п����������� 
    �����뵼��finish ָ����� end_ofstorage ʱ������Ҫ���� 
    */  
      
    iterator i = start;//old_start  
    iterator new_finish;//�����ݵĻ�ָ��old_start�����ݵĻ�ָ��new_start  
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
      
    /*ԭʼ���ݵİ���*/  
    while(i<position) { *new_finish++ = *i++; }  
    /*���벿�ֵĳ�ʼ��*/  
    while(n--) { *new_finish++ = x; }  
    /*ԭʼ���ݵİ���*/  
    while(i<finish) { *new_finish++ = *i++; }  
    finish = new_finish;  
  
    /*������Ҫ�ͷ�ԭ�пռ�*/  
    if(needToDestory==true){  
        while(old_start != old_finish){  
            old_start->~T();  
            old_start++;  
        }  
    }  
}  
#endif

#endif /* #ifndef __my_vector_h__ */  
