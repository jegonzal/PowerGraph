#!/usr/bin/env ruby

require 'json'

class String
  alias :old_to_f :to_f
  def to_f
    return 1e99 if 0 >= self.length
    self.old_to_f
  end
end

class Handler

  def gather_edges(json)
    {:edges => :IN_EDGES}
  end

  def scatter_edges(json)
    {:edges => :OUT_EDGES}
  end

  def gather(json)
    edge = json["params"]["edge"]
    # get distance to this vertex thru edge
    edge_dist = edge["state"].to_f
    nbr_dist = edge["source"]["state"].to_f
    vertex_dist = edge_dist + nbr_dist
    # return gather_type
    {:result => vertex_dist.to_s}
  end

  def merge(json)
    # take min of two gathers
    left = json["state"].to_f
    right = json["params"]["other"].to_f
    {:result => [left, right].min.to_s}
  end

  def apply(json)
    # take min of vertex dist and gather dist
    vertex_dist = json["params"]["vertex"]["state"].to_f
    gather_dist = json["params"]["gather"].to_f
    return {:vertex => gather_dist.to_s} if (gather_dist < vertex_dist)
    return {}
  end

  def scatter(json)
    edge = json["params"]["edge"]
    # get distance to neighbor thru this edge
    vertex_dist = json["params"]["vertex"]["state"].to_f
    edge_dist = edge["state"].to_f
    nbr_dist = edge["target"]["state"].to_f
    return {:signal => :TARGET} if (nbr_dist > (vertex_dist + edge_dist))
    return {}
  end

end

h = Handler.new

begin

  # first line tells us some many bytes
  bytes = STDIN.readline.to_i

  json = JSON.parse(STDIN.read bytes)
  method = json["method"]

  # invoke
  raise IOError, "Missing method field" if (nil == method)
  break if "exit" == method
  raise IOError, "Unrecognized method: " + method if not h.respond_to? method
  json = h.send method, json

  # return
  return_str = JSON.generate json
  puts return_str.length
  print return_str
  STDOUT.flush

rescue EOFError
  STDERR.write "Pipe broken"
end while true
