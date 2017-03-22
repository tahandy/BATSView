function o_ui32 = morton2(x_ui32, y_ui32)

    o_ui32 = Part1By1(x_ui32) + bitsll(Part1By1(y_ui32),1);