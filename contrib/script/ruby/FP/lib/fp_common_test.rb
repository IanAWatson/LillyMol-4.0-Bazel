require_relative "fp_common"
require "test/unit"

class TestParseFpToken < Test::Unit::TestCase

  def test_empty
    fp = ''
    legacy, radius, atype = FpCommon.parse_fp_token(fp)
    assert_nil(legacy)
    assert_nil(radius)
    assert_nil(atype)
  end

  [['1', 'C', ['C', nil, nil]],
   ['2', 'SB', ['SB', nil, nil]],
   ['3', 'ABC', ['ABC', nil, nil]],
   ['4', '1', [nil, '1', nil]],
   ['5', '2', [nil, '2', nil]],
   ['6', '10', [nil, '10', nil]],
   ['7', 'C1', ['C', '1', nil]],
   ['8', 'C10', ['C', '10', nil]],
   ['9', 'ABC10', ['ABC', '10', nil]],
   ['10', ':xyz', [nil, nil, 'xyz']],
   ['11', '3:xyz', [nil, '3', 'xyz']],
   ['12', '23:xyz', [nil, '23', 'xyz']],
   ['13', 'Y23:xyz', ['Y', '23', 'xyz']],
   ['13', 'YY23:xyz', ['YY', '23', 'xyz']]
  ].each do |scenario, initial, expected|
    test_name = "test_#{scenario}"
    define_method(test_name) {
      result = FpCommon.parse_fp_token(initial)
      assert_equal(result, expected)
    }
  end
end
