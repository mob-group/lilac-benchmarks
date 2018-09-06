#include <unistd.h>
#include <signal.h>
#include <sys/mman.h>
#include <vector>

class ReadLocation
{
    static struct sigaction old_sigaction;

    static void handler(int sig, siginfo_t *si, void *unused)
    {
        bool any_hit = false;

        for(auto* location : trapped_locations)
        {
            any_hit |= location->hit(si->si_addr);
        }

        if(any_hit) return;

        old_sigaction.sa_sigaction(sig, si, unused);
    }

public:
    ReadLocation(ReadLocation&) = delete;

    ReadLocation() : begin(nullptr), size(0), aligned_begin(0),
                 aligned_end(0), changed(true)
    {
        if(trapped_locations.empty())
        {
            struct sigaction action;
            action.sa_sigaction = handler;
            sigemptyset(&action.sa_mask);
            action.sa_flags = SA_SIGINFO;
            sigaction(SIGSEGV, &action, &old_sigaction);
        }

        trapped_locations.push_back(this);
    }

    void set(char* p, size_t n)
    {
        if(p != begin || n != size)
        {
            exit();
            
            size_t pagesize = getpagesize();

            begin = aligned_begin = (char*)p;
            aligned_end = p+n;
            size  = n;

            padding = {};
            rightpadding = {};

            while(((long)aligned_begin) % pagesize)
            {
                if(aligned_begin < aligned_end)
                    padding.push_back(*aligned_begin);
                aligned_begin ++;
            }

            while(((long)aligned_end+1) % pagesize)
            {
                if(aligned_end > aligned_begin)
                    rightpadding.push_back(*(aligned_end-1));
                aligned_end --;
            }

            init();
        }

        update_if_necessary();
    }

    virtual void init() { };
    virtual void update() { };
    virtual void exit() { };

    void update_if_necessary()
    {
        if(has_changed())
        {
            update();

            changed = false;

            for(size_t i = 0; i < padding.size(); i++)
                padding[i] = begin[i];

            for(size_t i = 0; i < rightpadding.size(); i++)
                rightpadding[i] = begin[size-i-1];
            
            if(aligned_begin < aligned_end)
                mprotect(aligned_begin, aligned_end-aligned_begin,
                         PROT_READ | PROT_EXEC);
        }
    }


    void clear_traps()
    {
        changed = true;
        if(aligned_begin < aligned_end)
            mprotect(aligned_begin, aligned_end-aligned_begin,
                     PROT_READ | PROT_WRITE | PROT_EXEC);
    }

    bool hit(void* ptr)
    {
        if(ptr >= aligned_begin && ptr < aligned_end)
        {
            if(!changed)
            {
                clear_traps();
            }
            return true;
        }

        return false;
    }

    bool has_changed()
    {
        if(changed) return true;

        for(size_t i = 0; i < padding.size(); i++)
            if(padding[i] != begin[i])
                return true;

        for(size_t i = 0; i < rightpadding.size(); i++)
            if(rightpadding[i] != begin[size-i-1])
                return true;

        return false;
    }


    static std::vector<ReadLocation*> trapped_locations;

protected:
    char*  begin;
    size_t size;

private:
    char*  aligned_begin;
    char*  aligned_end;
    bool   changed;

    std::vector<char> padding;
    std::vector<char> rightpadding;
};

std::vector<ReadLocation*> ReadLocation::trapped_locations{};
struct sigaction ReadLocation::old_sigaction{};


template<typename Tin, typename Tout,
         void(*u)(Tin*,int,Tout&),
         void(*construct)(int,Tout&) = nullptr,
         void(*destruct)(int,Tout&) = nullptr>
class ReadObject : public ReadLocation
{
public:
    const Tout& operator()(Tin* p, size_t n)
    {
        begin = p;
        size = n;
        ReadLocation::set((char*)p, sizeof(Tin)*n);
        return value;
    }

    void init() final
    {
        if(construct) construct(this->size, this->value);
    }

    void update() final
    {
        if(u) u(this->begin, this->size, this->value);
    }

    void exit() final
    {
        if(destruct) destruct(this->size, this->value);
    }

protected:
    Tin*   begin;
    size_t size;
    Tout   value;
};

template<typename Tin,typename Tout,
         void(*update)(Tin*,int,Tout&),
         void(*construct)(int,Tout&) = nullptr,
         void(*destruct)(int,Tout&) = nullptr>
class WriteObject
{
private:
    class Scope
    {
    public:
        Scope(Tin* b, Tout& p, size_t n) : begin(b), dev_ptr(p), size(n) { }

        ~Scope()
        {
            update(begin, size, dev_ptr);
        }

        operator Tout ()
        {
            return dev_ptr;
        }

    private:

        Tin* begin;
        Tout& dev_ptr;
        size_t size;
    };
public:
    WriteObject() : size(0) { }

    Scope operator()(Tin* p, size_t n)
    {
        if(n != size)
        {
            if(size != 0) destruct(size, dev_ptr);
            size = n;
            if(size != 0) construct(size, dev_ptr);
        }

        return Scope{p, dev_ptr, n};
    }

private:
    Tout dev_ptr;
    size_t size;
};
