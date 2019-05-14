#ifndef _Int_Stack_h__
#define _Int_Stack_h__

class IntStack
{
public:
    IntStack();
	IntStack(int maxStackSize);
    ~IntStack();

    void push(int elem);
    void pop();
    int top();
    bool empty();
	void reset();
private:
    int *elems_;
    int top_;
	int maxSize_;
};
#endif /* _Int_Stack_h__ */