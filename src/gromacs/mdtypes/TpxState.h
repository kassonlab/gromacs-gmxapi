//
// Created by Eric Irrgang on 8/27/17.
//

#ifndef GMX_MDTYPES_TPXSTATE_H
#define GMX_MDTYPES_TPXSTATE_H

#include <memory>

struct t_inputrec;
class t_state;
struct gmx_mtop_t;

namespace gmx
{

class TpxState
{
    private:
        std::shared_ptr<t_inputrec>      inputrecInstance_;
        std::shared_ptr<t_state>         stateInstance_;
        gmx_mtop_t               *mtop_;
        bool                      initialized_;
        bool                      dirty_;
    public:
        TpxState();
        ~TpxState();

        // Copy semantics TBD
        TpxState(const TpxState&) = delete;
        TpxState& operator=(const TpxState&) = delete;

        // Move should be okay
        TpxState(TpxState&&) noexcept = default;
        TpxState& operator=(TpxState&&) noexcept = default;

        static std::unique_ptr<TpxState> initializeFromFile(const char* filename);
        static std::unique_ptr<TpxState> initializeFromFile(const std::string& filename);

        t_inputrec* getRawInputrec();
        gmx_mtop_t* getRawMtop();
        t_state*    getRawState();

        bool isInitialized() const;

};


} // end namespace gmx

#endif //GROMACS_TPXSTATE_H
