
namespace kg
{
namespace io
{

inline IOAdios2::IOAdios2() : ad_{MPI_COMM_WORLD, adios2::DebugON} {}

inline File IOAdios2::openFile(const std::string& name, const Mode mode,
                               MPI_Comm comm, const std::string& io_name)
{
  return File{new FileAdios2{ad_, name, mode, io_name}};
}

inline Engine IOAdios2::open(const std::string& name, const Mode mode,
                             MPI_Comm comm, const std::string& io_name)
{
  return {openFile(name, mode, comm, io_name), comm};
}

} // namespace io
} // namespace kg
