#!/usr/bin/env ruby

require 'json'

RESET_PROB = 0.15;

class Handler

  def gather_edges(json)
    {:edges => :IN_EDGES}
  end

  def scatter_edges(json)
    last_change = json['state'].to_f
    return {:edges => :OUT_EDGES} if last_change > 1e-2
    return {:edges => :NO_EDGES}
  end

  def gather(json)
    edge = json['params']['edge']
    edge_weight = edge['state'].to_f            # normalized edge weight
    nbr_rank = edge['source']['state'].to_f
    {:result => (edge_weight * nbr_rank).to_s}
  end

  def merge(json)
    # take sum of two gathers
    left = json['state'].to_f
    right = json['params']['other'].to_f
    {:result => (left + right).to_s}
  end

  def apply(json)
    newval = json['params']['gather'].to_f + RESET_PROB
    last_change = (newval - json['params']['vertex']['state'].to_f).abs
    return {:vertex => newval.to_s, :program => last_change.to_s}
  end

  def scatter(json)
    last_change = json['state'].to_f
    return {:signal => :TARGET} if last_change > 1e-2
    return {}
  end
  
  def init_edge(json)
    source = json['params']['edge']['source']
    weight = (1.0 - RESET_PROB) / source['num_out_edges'].to_f
    return {:edge => weight.to_s}
  end
  
  def init_vertex(json)
    return {:vertex => 1.0.to_s}
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
  
  json = if h.respond_to? method
    h.send method, json
  else
    {}
  end

  # return
  return_str = JSON.generate json
  puts return_str.length
  print return_str
  STDOUT.flush

rescue EOFError
  STDERR.write 'Pipe broken'
  exit
end while true
