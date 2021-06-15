

#include <string>
#include <list>
namespace outputWriter{

    class pvdWriter{

    public:pvdWriter();
       virtual ~pvdWriter();
        void writeFile(const std::string &filename, std::list<std::string> str );

    };
}