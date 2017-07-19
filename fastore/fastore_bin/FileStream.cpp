/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#if defined (_WIN32)
#	define _CRT_SECURE_NO_WARNINGS
#	pragma warning(disable : 4996) // D_SCL_SECURE
#	pragma warning(disable : 4244) // conversion uint64 to uint32
//#	pragma warning(disable : 4267)
#	define FOPEN	fopen
#	define FSEEK	_fseeki64
#	define FTELL	_ftelli64
#	define FCLOSE	fclose
#elif __APPLE__	// Apple by default suport 64 bit file operations (Darwin 10.5+)
#	define FOPEN	fopen
#	define FSEEK	fseek
#	define FTELL	ftell
#	define FCLOSE	fclose
#else
#	if !defined(_LARGEFILE_SOURCE)
#		define _LARGEFILE_SOURCE
#		if !defined(_LARGEFILE64_SOURCE)
#			define _LARGEFILE64_SOURCE
#		endif
#	endif
#	if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#		undef _FILE_OFFSET_BITS
#	endif
#	if !defined(_FILE_OFFSET_BITS)
#		define _FILE_OFFSET_BITS 64
#	endif
#	define FOPEN	fopen64
#	define FSEEK	fseeko64
#	define FTELL	ftello64
#	define FCLOSE	fclose
#endif

#include "Globals.h"

#include <stdio.h>

// for memory mapped files -- are we using this ????
//
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <errno.h>
#include <string.h>
//
//

#include "FileStream.h"
#include "Exception.h"


struct IFileStream::FileStreamImpl
{
	FILE* file;

	FileStreamImpl()
		:	file(NULL)
	{}
};


IFileStream::IFileStream()
{
	impl = new FileStreamImpl();
}


IFileStream::~IFileStream()
{
	ASSERT(impl != NULL);
	delete impl;
}


void IFileStream::SetBuffering(bool enable_)
{
	if (enable_)
		setvbuf(impl->file, NULL, _IOFBF, 1 << 23);
	else
		setvbuf(impl->file, NULL, _IONBF, 0);
}


struct IMemoryStream::MemoryStreamImpl
{
	int32 fileDescriptor;
	void* memory;
	uint64 size;

	MemoryStreamImpl()
		:	fileDescriptor(0)
		,	memory(NULL)
		,	size(0)
	{}
};


IMemoryStream::IMemoryStream()
{
	impl = new MemoryStreamImpl();
}


IMemoryStream::~IMemoryStream()
{
	delete impl;
}


void* IMemoryStream::Pointer()
{
	return impl->memory;
}


FileStreamReader::FileStreamReader(const std::string& fileName_)
	:	size(0)
	,	position(0)
{
	FILE* f = FOPEN(fileName_.c_str(), "rb");
	if (f == NULL)
	{
		throw Exception("Cannot open file to read: " + fileName_);
	}
	impl->file = f;

	FSEEK(impl->file, 0, SEEK_END);
	size = FTELL(impl->file);

	FSEEK(impl->file, 0, SEEK_SET);
	position = 0;
}


FileStreamReader::~FileStreamReader()
{
	ASSERT(impl != NULL);

	if (impl->file != NULL)
		FCLOSE(impl->file);
}


void FileStreamReader::Close()
{
	ASSERT(impl->file != NULL);

	FCLOSE(impl->file);
	impl->file = NULL;
}


int64 FileStreamReader::Read(uchar *mem_, uint64 size_)
{
	int64 n = fread(mem_, 1, size_, impl->file);
	if (n >= 0)
		position += n;
	return n;
}


void FileStreamReader::SetPosition(uint64 pos_)
{
	ASSERT(impl->file != NULL);

	if (pos_ == position)
		return;

	if (pos_ > size)
		throw Exception("Position exceeds stream size");

	FSEEK(impl->file, pos_, SEEK_SET);
	position = pos_;
}


MemoryStreamReader::MemoryStreamReader(const std::string &fileName_)
	:	position(0)
	,	memoryBuffer(NULL)
{
	int32 fd = open(fileName_.c_str(), O_RDONLY);

	if (fd < 0)
		throw Exception("Cannot open file: " + fileName_);

	struct stat s;
	if (fstat(fd, &s) < 0 || s.st_size == 0)
	{
		close(fd);
		throw Exception("Cannot stat file or file is empty.");
	}

	const void* region = mmap(NULL, s.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (region == MAP_FAILED)
	{
		close(fd);
		throw Exception("Cannot mmap file: " + std::string(strerror(errno)));
	}

	impl->fileDescriptor = fd;
	impl->size = s.st_size;
	impl->memory = (void*)region;
	position = 0;

	memoryBuffer = new IBuffer((byte*)impl->memory, impl->size);
}


MemoryStreamReader::~MemoryStreamReader()
{
	Close();

	if (memoryBuffer != NULL)
		delete memoryBuffer;
}


void MemoryStreamReader::Close()
{
	if (impl->memory != NULL)
	{
		munmap(impl->memory, impl->size);
		impl->memory = NULL;
	}

	if (impl->fileDescriptor > 0)
	{
		close(impl->fileDescriptor);
		impl->fileDescriptor = 0;
	}

	if (memoryBuffer != NULL)
	{
		delete memoryBuffer;
		memoryBuffer = NULL;
	}

	impl->size = 0;
	position = 0;
}


void MemoryStreamReader::SetPosition(uint64 pos_)
{
	if (pos_ > impl->size)
		throw Exception("Position exceeds stream size");

	position = pos_;
}


uint64 MemoryStreamReader::Size() const
{
	return impl->size;
}


int64 MemoryStreamReader::Read(uchar *mem_, uint64 size_)
{
	int64 toRead = size_;
	if (position + toRead > impl->size)
		toRead = impl->size - position;

	if (toRead > 0)
	{
		std::copy((const uchar*)impl->memory + position, (const uchar*)impl->memory + position + toRead, mem_);
		position += toRead;
	}
	return toRead;
}


int64 MemoryStreamReader::Attach(uchar *&mem_, uint64 size_)
{
	int64 toRead = size_;
	if (position + toRead > impl->size)
		toRead = impl->size - position;

	if (toRead > 0)
	{
		mem_ = (uchar*)impl->memory + position;
		position += toRead;
	}
	return toRead;
}


FileStreamWriter::FileStreamWriter(const std::string& fileName_)
	:	position(0)
{
	FILE* f = FOPEN(fileName_.c_str(), "wb");
	if (f == NULL)
	{
		throw Exception("Cannot open file to write: " + fileName_);
	}
	impl->file = f;
}


FileStreamWriter::~FileStreamWriter()
{
	ASSERT(impl != NULL);

	if (impl->file != NULL)
		FCLOSE(impl->file);
}


void FileStreamWriter::Close()
{
	ASSERT(impl->file != NULL);

	FCLOSE(impl->file);
	impl->file = NULL;
}


int64 FileStreamWriter::Write(const uchar *mem_, uint64 size_)
{
	int64 n = fwrite(mem_, 1, size_, impl->file);
	if (n >= 0)
		position += n;
	return n;
}


void FileStreamWriter::SetPosition(uint64 pos_)
{
	ASSERT(impl->file != NULL);

	FSEEK(impl->file, pos_, SEEK_SET);
	position = pos_;
}

IMultiFileStreamReader::IMultiFileStreamReader(const std::vector<std::string> &fileNames_, const IFileFuncImpl* fileFuncImpl_)
	:	fileNames(fileNames_)
	,	fileFuncImpl(fileFuncImpl_)
	,	idx(0)
{
	ASSERT(fileFuncImpl_ != NULL);

	if (fileNames_.size() == 0)
		throw Exception("Empty file list.");

	OpenNext();
}


IMultiFileStreamReader::~IMultiFileStreamReader()
{
	Close();

	if (fileFuncImpl != NULL)
		delete fileFuncImpl;
}


void IMultiFileStreamReader::Close()
{
	if (impl->file != NULL)
	{
		fileFuncImpl->Close(impl->file);
		impl->file = NULL;
	}
}


bool IMultiFileStreamReader::OpenNext()
{
	if (impl->file != NULL)
	{
		fileFuncImpl->Close(impl->file);
		impl->file = NULL;
	}

	if (idx >= fileNames.size())
		return false;

	FILE* f = (FILE*)fileFuncImpl->Open(fileNames[idx].c_str(), "rb");
	if (f == NULL)
		throw Exception(("Cannot open file to read:" + fileNames[idx]).c_str());

	idx++;
	impl->file = f;

	return true;
}


int64 IMultiFileStreamReader::Read(uchar *mem_, uint64 size_)
{
	if (impl->file == NULL)
		return 0;

	int64 n = fileFuncImpl->Read(impl->file, mem_, size_);

	while (n < (int64)size_ && OpenNext())
	{
		uint64 n2 = fileFuncImpl->Read(impl->file, mem_ + n, size_ - n);
		if (n2 > 0)
			n += n2;
	}

	return n;
}


struct StdFileFuncImpl : public IFileFuncImpl
{
	void* Open(const char* filename_, const char* flags_) const
	{
		FILE* f = FOPEN(filename_, flags_);
		return f;
	}

	void Close(void* file_) const
	{
		FCLOSE((FILE*)file_);
	}

	int64 Read(void* file_, byte* mem_, uint64 size_) const
	{
		return fread(mem_, 1, size_, (FILE*)file_);
	}
};


#include <zlib.h>

struct GzFileFuncImpl : public IFileFuncImpl
{
	void* Open(const char* filename_, const char* flags_) const
	{
		gzFile f = gzopen(filename_, flags_);
		if (f != NULL)
		{
			gzbuffer(f, 32 << 20);		// 32 MB for the buffer, internally it's 3x ; TODO: parameter
		}
		return f;
	}

	void Close(void* file_) const
	{
		gzclose((gzFile)file_);
	}

	int64 Read(void* file_, byte* mem_, uint64 size_) const
	{
		return gzread((gzFile)file_, (voidp)mem_, size_);
	}
};


MultiFileStreamReader::MultiFileStreamReader(const std::vector<std::string> &fileNames_)
	: IMultiFileStreamReader(fileNames_, new StdFileFuncImpl())
{}


MultiFileStreamReaderGz::MultiFileStreamReaderGz(const std::vector<std::string> &fileNames_)
	: IMultiFileStreamReader(fileNames_, new GzFileFuncImpl())
{}

