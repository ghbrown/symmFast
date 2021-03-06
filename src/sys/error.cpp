#include <sys/sys.hpp>

#include <petscsys.h>

#include <tuple>
#include <queue>
#include <string>
#include <iostream>

namespace sf
{

namespace detail
{

#define SF_ERROR_INITIAL -1
#define SF_ERROR_REPEAT  -2

static constexpr const char error_bar[] = "-----------------------------------------------------------------------------------";

class backtrace
{
private:
  using entry_type = std::tuple<std::string,std::string,int>;
  using stack_type = std::queue<entry_type>;

  const sf_error_t error_;
  stack_type       stack_;

public:
  backtrace(sf_error_t errc) noexcept : error_(errc) { }

  const entry_type& last() const noexcept { return stack_.back(); }

  sf_error_t error_code() const noexcept  { return error_;        }

  void push(std::string file, std::string func, int line) noexcept
  {
    stack_.push(std::make_tuple(std::move(file),std::move(func),line));
  }

  void unwind() noexcept
  {
    int i = 0;
    std::cout<<"SymmFast error "<<error_bar<<std::endl;
    std::cout<<"Error code "<<error_<<std::endl;
    while (!stack_.empty()) {
      int         line;
      std::string func,file;

      std::tie(file,func,line) = stack_.front();
      std::cout<<"ERROR: ["<<i++<<"] "<<func<<"() line "<<line<<" in file "<<file<<std::endl;
      stack_.pop();
    }
  }
};

sf_error_t error_handler(const char *file, const char *func, int line, MPI_Comm comm, int errc)
{
  static backtrace bt(errc);

  bt.push(file,func,line);
  if (std::get<1>(bt.last()) == "main") {
    // char errorstr[MPI_MAX_ERROR_STRING];
    // int  resultlen; // unused

    // MPI_Error_string(errc,static_cast<char*>(errorstr),&resultlen);
    bt.unwind();
    MPI_Abort(comm,bt.error_code()); // noreturn
    __builtin_unreachable();
  }
  return SF_ERROR_REPEAT;
}

} // namespace detail

} // namespace sf
