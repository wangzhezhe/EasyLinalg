#include <iostream>
#include <assert.h>

template <typename T, uint Size>
class Vec
{
public:
    static const int NUM_COMPONENTS = Size;
    T Components[NUM_COMPONENTS];
    Vec() {}
    Vec(T init)
    {
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            this->Components[i] = init;
        }
    }

    Vec(const Vec &src)
    {
        assert(this->NUM_COMPONENTS == src.NUM_COMPONENTS);
        for (uint i = 0; i < Size; ++i)
        {
            this->Components[i] = src[i];
        }
    }

    Vec operator+(const Vec &t) const
    {
        assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
        Vec<T, Size> result;
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            result[i] = this->Components[i] + t.Components[i];
        }
        return result;
    }

    Vec operator-(const Vec &t) const
    {
        assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
        Vec<T, Size> result;
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            result[i] = this->Components[i] - t.Components[i];
        }
        return result;
    }

    T &operator[](int index)
    {
        assert(index >= 0);
        assert(index < this->NUM_COMPONENTS);
        return this->Components[index];
    }

    const T &operator[](int index) const
    {
        assert(index >= 0);
        assert(index < this->NUM_COMPONENTS);
        return this->Components[index];
    }

    Vec &operator=(const Vec &t)
    {
        assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            this->Components[i] = t.Components[i];
        }
        return *this;
    }
    void Show()
    {
        for (int i = 0; i < this->NUM_COMPONENTS; i++)
        {
            printf("%lf ", 1.0 * this->Components[i]);
        }
        printf("\n");
    }
};

int main()
{
    Vec<double, 3> v1;
    for (int i = 0; i < 3; i++)
    {
        v1[i] = 2.0;
    }
    Vec<double, 3> v2;
    for (int i = 0; i < 3; i++)
    {
        v2[i] = 3.0;
    }
    auto vadd = v1 + v2;
    std::cout << "test";
    vadd.Show();
    return 0;
}