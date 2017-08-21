#include "gmxapi/md.h"
#include <gtest/gtest.h>

#include "gmxapi/runner.h"
#include "md-impl.h"
#include "gromacs/utility/keyvaluetree.h"

TEST(ApiModuleMD, Construction)
{
    {   // Check default construction and destruction
        gmxapi::MDProxy proxy{};
    }
    auto module = gmxapi::mdFromTpr("topol.tpr");
    // Helper function should return a non-null unique_ptr
    ASSERT_TRUE(module.get() != nullptr);
    ASSERT_EQ(module->info(), "MDStatePlaceholder initialized with filename: \"topol.tpr\"\n");
}

TEST(ApiModuleMD, Binding)
{
    class MyRunner : public gmxapi::IMDRunner
    {
        private:
            std::unique_ptr<gmxapi::MDBuilder> builder_;
        public:
            void registerMDBuilder(std::unique_ptr<gmxapi::MDBuilder> builder) override
            {
                builder_ = std::move(builder);
            };
            gmxapi::Status run() override
            {
                if (builder_ != nullptr)
                {
                    builder_->build();
                }
//                else assert(false);
                return gmxapi::Status();
            };
    };

    {
        // Dummy state object to test registration protocol.
        class MyDummyMD: public gmxapi::MDState
        {
            public:
                MyDummyMD() = default;
                virtual ~MyDummyMD() = default;

                std::unique_ptr<gmxapi::MDBuilder> builder() override
                {
                    // Todo: rewrite!!!
                    return nullptr;
                }
        };
        auto runner = MyRunner();
        MyDummyMD module{};
        runner.registerMDBuilder(module.builder());
        runner.run();
    }


    // Register MD with a runner
    auto module = gmxapi::mdFromTpr("topol.tpr");
    auto runner = MyRunner();
    module->bind(&runner);
    runner.run();
}
