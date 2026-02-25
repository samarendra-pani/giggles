#ifndef TEST_COMPONENTFINDER_H
#define TEST_COMPONENTFINDER_H

#include "../componentfinder.h"
#include "../../../tests_data.h"
#include <cassert>

void test_componentfinder() {

    const std::vector<uint32_t> values = {100, 200, 300, 400, 500};
    ComponentFinder<uint32_t> componentfinder(values);

    componentfinder.merge(100, 300);

    assert(componentfinder.find(100) == 100);
    assert(componentfinder.find(200) == 200);
    assert(componentfinder.find(300) == 100);
    assert(componentfinder.find(400) == 400);
    assert(componentfinder.find(500) == 500);

    componentfinder.merge(200, 400);

    assert(componentfinder.find(100) == 100);
    assert(componentfinder.find(200) == 200);
    assert(componentfinder.find(300) == 100);
    assert(componentfinder.find(400) == 200);
    assert(componentfinder.find(500) == 500);

    componentfinder.merge(300, 500);

    assert(componentfinder.find(100) == 100);
    assert(componentfinder.find(200) == 200);
    assert(componentfinder.find(300) == 100);
    assert(componentfinder.find(400) == 200);
    assert(componentfinder.find(500) == 100);

    componentfinder.merge(200, 500);

    assert(componentfinder.find(100) == 100);
    assert(componentfinder.find(200) == 100);
    assert(componentfinder.find(300) == 100);
    assert(componentfinder.find(400) == 100);
    assert(componentfinder.find(500) == 100);

    assert_msg(true, "ComponentFinder", "Component creation with 5 values");
}


#endif // TEST_COMPONENTFINDER_H