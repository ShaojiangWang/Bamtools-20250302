// ***************************************************************************
// HostInfo_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides DNS lookup functionality for hostname/IP addresses
// ***************************************************************************

#ifndef HOSTINFO_P_H
#define HOSTINFO_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/internal/io/HostAddress_p.h"
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

class HostInfo {

    public:
        enum ErrorType { NoError = 0
                       , HostNotFound
                       , UnknownError
                       };

    // ctors & dtor
    public:
        HostInfo(void);
        HostInfo(const HostInfo& other)
	    : m_hostName(other.m_hostName), m_addresses(other.m_addresses), 
              m_error(other.m_error), m_errorString(other.m_errorString)
	{ }
        HostInfo(HostInfo&& other) 
	    : m_hostName(std::move(other.m_hostName)), 
              m_addresses(std::move(other.m_addresses))
	      ,m_error(other.m_error), 
              m_errorString(std::move(other.m_errorString))
	{}
        HostInfo& operator=(const HostInfo& other);
        HostInfo& operator=(HostInfo&& other);
        ~HostInfo(void);

    // HostInfo interface
    public:
        std::string HostName(void) const;
        void SetHostName(const std::string& name);

        std::vector<HostAddress> Addresses(void) const;
        void SetAddresses(const std::vector<HostAddress>& addresses);

        HostInfo::ErrorType GetError(void) const;
        std::string GetErrorString(void) const;

    // internal methods
    private:
        void SetError(const HostInfo::ErrorType error);
        void SetErrorString(const std::string& errorString);

    // static methods
    public:
        static HostInfo Lookup(const std::string& hostname,
                               const std::string& port);

    // data members
    private:
        std::string m_hostName;
        std::vector<HostAddress> m_addresses;
        HostInfo::ErrorType m_error;
        std::string m_errorString;
};

} // namespace Internal
} // namespace BamTools

#endif // HOSTINFO_P_H
