%%defining the source-receiver class
classdef SourceReceiver
    properties 
        ReceiverCoord = [0,0,0];
        SourceCoord;
    end
    methods
        function obj = SourceReceiver(SCoord, RCoord)
            obj.ReceiverCoord = RCoord;
            obj.SourceCoord = SCoord;
        end
        function f = get_source()
            f = obj.SourceCoord;
        end
        function f = get_receiver()
            f=obj.ReceiverCoord;
        end
    end
end

