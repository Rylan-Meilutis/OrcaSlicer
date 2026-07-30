#pragma once

#include <string>
#include <vector>

namespace Slic3r {

class AppConfig;

struct ProfileSource
{
    enum class Format { Orca, Prusa };

    std::string id;
    std::string name;
    std::string url;
    Format      format {Format::Prusa};
    long long   last_sync {0};
    bool        enabled {false};
};

struct ProfileSourceSyncResult
{
    size_t      printers {0};
    size_t      filaments {0};
    size_t      processes {0};
    size_t      skipped_options {0};
    std::string error;

    bool success() const { return error.empty(); }
};

// Stores source definitions in AppConfig and installs each source into its own
// local preset bundle. Imported preset names are qualified by source (and by
// vendor for native repositories), preventing name-based preset lookup from
// confusing otherwise identical profiles from different upstreams. This
// ownership boundary also makes source removal safe: user presets and presets
// installed by another source are never touched.
class ProfileSourceManager
{
public:
    explicit ProfileSourceManager(AppConfig &config);

    std::vector<ProfileSource> sources() const;
    void                       set_sources(const std::vector<ProfileSource> &sources);
    void                       add(const ProfileSource &source);
    bool                       remove(const std::string &id, std::string &error);
    ProfileSourceSyncResult    sync(const ProfileSource &source);
    std::vector<ProfileSource> stale_enabled_sources(long long max_age_seconds = 24 * 60 * 60) const;

    static std::vector<ProfileSource> prusa_sources();
    static std::vector<ProfileSource> default_sources();
    static std::string                make_id(const std::string &name, const std::string &url);
    // Converts an extracted PrusaSlicer/SuperSlicer profile tree. Public so
    // conversion can be regression-tested without network access.
    static ProfileSourceSyncResult    convert_prusa_profiles(const std::string &input_root,
                                                              const std::string &output_root,
                                                              const std::string &source_namespace = {});
    static ProfileSourceSyncResult    convert_orca_profiles(const std::string &input_root,
                                                             const std::string &output_root,
                                                             const std::string &source_id,
                                                             const std::string &source_name);

private:
    AppConfig &m_config;
};

} // namespace Slic3r
