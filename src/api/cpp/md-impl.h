#ifndef GMXAPI_MD_IMPL_H
#define GMXAPI_MD_IMPL_H
/*! \file
 * \brief Declarations for molecular dynamics API implementation details.
 *
 * \ingroup gmxapi
 */

#include "gmxapi/md.h"
#include <memory>
#include <string>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

// Container for MD engine input parameters.
class MDInput
{
    public:
        MDInput();

        /// Take over and wrap input data structures.
        MDInput(std::unique_ptr<t_inputrec> &&inputRecord,
                std::unique_ptr<t_state>    &&state,
                std::unique_ptr<gmx_mtop_t> &&topology);

        /// Get input record and state from TPR file.
        static std::unique_ptr<MDInput> from_tpr_file(std::string filename);

        int nAtoms() const;

        /// Return a copy of the KeyValueTreeObject
        gmx::KeyValueTreeObject params() const;

        const t_state* state() const;

        std::unique_ptr<t_inputrec> inputRecord_;
        std::unique_ptr<t_state>    state_;
        std::unique_ptr<gmx_mtop_t> topology_;
};

/*!
 * \brief Base class for the run-time objects fronted by MDProxy.
 *
 * An MD task can have a handle before, during, or after execution, and the local handle may refer
 * to a different implementation class depending on whether execution takes place locally or remotely,
 * and whether distributed data structures are cached locally, etc.
 *
 * Pure virtual base class is not instantiated by clients directly. Objects are created
 * by other API objects or helper functions. \see mdFromTpr()
 *
 * A State instance provides the instantaneous implementation for the owning object.
 *
 * The MDState implementation is responsible for implementing the bind() method, allowing a
 * runner proxy and MD proxy to be translated into an actual runner and MD Engine instance.
 *
 * \todo For ModuleMD to be user-extensible, there needs to be a public interface for
 * the State classes as well as the proxy classes.
 */
class MDState
{
    public:
        MDState() = default;
        virtual ~MDState();
        MDState(const MDState&) = delete;
        MDState& operator=(const MDState&) = delete;
        /*!
         * \brief Get a builder for an MD Engine from this proxy object
         *
         * \return ownership of a MD Engine builder implementing the gmxapi::ModuleBuilder interface.
         */
        virtual std::unique_ptr<MDBuilder> builder() = 0;

        /// Allow implementing classes to provide information in a generic way.
        virtual std::string info();
};


/*! \brief Data based on C structs.
 *
 * Provide a ModuleMD implementation that can be created from a MDInput struct
 * of the sort produced by MDInput::from_tpr_file().
 */
class MDStateFromMDInput : public MDState
{
    private:
        std::unique_ptr<MDInput> input_;
        std::string metadata_;

    public:
        MDStateFromMDInput();
        virtual ~MDStateFromMDInput() final;

        explicit MDStateFromMDInput(std::unique_ptr<MDInput> input);
        MDStateFromMDInput(std::unique_ptr<MDInput> input, std::string metadata);
        virtual std::unique_ptr<MDBuilder> builder() override;
        virtual std::string info() override;
};

/*!
 * \brief Thin implementation just for holding a tpr filename
 */
class MDStatePlaceholder : public MDState
{
    public:
        std::string filename_;
        explicit MDStatePlaceholder(const std::string& filename) : filename_{filename}
        {};

        virtual std::string info() override
        {
            std::string output("MDStatePlaceholder initialized");
            output.append(" with filename: \"");
            output.append(filename_);
            output.append("\"\n");
            return output;
        }

        std::unique_ptr<MDBuilder> builder() override
        {
            class NonBuilder : public MDBuilder
            {
                public:
                    std::string filename_;
                    explicit NonBuilder(const std::string& filename) : filename_{filename} {};
                    std::unique_ptr<MDEngine> build() override
                    {
                        return nullptr;
                    }

                    std::string inputAsTprFilename() override
                    {
                        return filename_;
                    }
            };
            std::unique_ptr<MDBuilder> builder = gmx::compat::make_unique<NonBuilder>(filename_);
            return builder;
        }
};

} // namespace gmxapi

#endif // header guard
